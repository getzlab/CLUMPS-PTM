# -- import packages: --------------------------------------------------------------------
import os
import contextlib
import tempfile
import subprocess
import glob
import pandas as pd
import numpy as np

AMINO_ACID_MAP = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E', \
         'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',\
         'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYX':'C',\
         'SEP':'S','TPO':'T','TYS':'Y','MK8':'L','M3L':'K','DPR':'P','DSN':'S',\
         'ALN':'A','DLY':'K','MVA':'V','MLE':'L','DLE':'L','CSO':'C','PTR':'Y',\
         'BMT':'T','DAL':'A','FAK':'K','MSE':'M','SMF':'A','HYP':'P'}

@contextlib.contextmanager
def gunzipper(gz_file: str):
    """
    Gunzipper: takes in a file name. Returns a valid file object in an unzipped form. Uses
    context manager to create a temporary file and save into that format. This
    was written for programmatic extraction of PDB files using Prody.

    Parameters:
    -----------
    gz_file [ required ]
        gzipped file
        type: str

    Returns:
    --------
    temporary file
        type: tempfile.NamedTemporaryFile

    """
    with tempfile.NamedTemporaryFile('r', suffix=os.path.splitext(gz_file)[1]) as temp_file:
        subprocess.check_call("gzip -dc {} >> {}".format(gz_file, temp_file.name), executable='/bin/bash', shell=True)
        yield temp_file

def add_corrected_fdr(results_df: pd.DataFrame, thresh_num: float = 0.1, weight_thresh_by_n: bool = False):
    """
    Creates corrected FDR for clumps-ptm results.
    Due to the sparsity of PTM modifications possible on a given
    protein, we perform FDR on the subset of structures that can
    yield a theoretical p-value < threshold

    Parameters:
    -----------
    results_df [ required ]
        results dataframe from clumps-ptm
        type: pd.DataFrame

    thresh_num [ required ]
        threshold for minimum theoretical p-value to consider for FDR
        type: float

    weight_thresh_by_n [ required ]
        whether to weight threshold for minimum theoretical p-value by hypothesis tested
        type: bool

    Returns:
    --------
    results_df with additional columns for corrected FDR
        type: pd.DataFrame

    """
    from scipy.special import comb
    from statsmodels.stats.multitest import multipletests

    corr_df = list()

    for sam in np.unique(results_df['clumpsptm_sampler']):
        _df_sam = results_df[results_df['clumpsptm_sampler']==sam].copy()

        _df_sam['fdr_max_pval'] = _df_sam.apply(lambda row: 1/comb(row['clumpsptm_sample_n'], row['clumpsptm_input_n']), 1)
        if weight_thresh_by_n:
            _df_sam['fdr_thresh'] = thresh_num / _df_sam.shape[0]
            _df_sam['fdr_pass'] = _df_sam['fdr_max_pval'] <= thresh_num / _df_sam.shape[0]
        else:
            _df_sam['fdr_thresh'] = thresh_num
            _df_sam['fdr_pass'] = _df_sam['fdr_max_pval'] <= thresh_num
        _df_sam['fdr_corr'] = 1

        _,_df_sam.loc[_df_sam['fdr_pass'],'fdr_corr'],_,_ = multipletests(
            _df_sam[_df_sam['fdr_pass']]['clumpsptm_pval'], method='fdr_bh',
        )

        corr_df.append(_df_sam)

    return pd.concat(corr_df)

def generate_clumpsptm_output(
    output_dir: str,
    protein_id: str = 'accession_number',
    site_id: str = 'variableSites',
    min_sites: int = 3,
    thresh_min_pval: bool = True,
    alphafold: bool = False
    ) -> pd.DataFrame:
    """
    Generate CLUMPS-PTM Output file.

    Parameters:
    -----------
    output_dir [ required ]
        output directory from clumps-ptm to pull results from
        type: str

    protein_id [ optional ]
        protein id used in input
        type: str

    site_id [ optional ]
        site id used in input
        type: str

    min_sites [ optional ]
        min-sites required for a given "clump"
        type: int

    thresh_min_pval [ optional ]
        whether to threshold empirical pvalues = 0 (due to insufficient permutations) to 1/n_permutations^e-1
        type: bool

    Returns:
    --------
    results dataframe
        type: pd.DataFrame

    """
    def _collapse_sites(df):
        """Collapser per feature type."""
        from statsmodels.stats.multitest import multipletests

        if alphafold:
            acc_df = df[
                ['geneSymbol','uniprot',protein_id,'clumpsptm_sampler',
                 'clumpsptm_niter','clumpsptm_pval','clumpsptm_input','clumpsptm_sample']
            ].set_index(protein_id).drop_duplicates().sort_values('clumpsptm_pval')
        else:
            acc_df = df[
                ['geneSymbol','uniprot',protein_id,'pdb','chain','clumpsptm_sampler',
                 'clumpsptm_niter','clumpsptm_pval','clumpsptm_input','clumpsptm_sample']
            ].set_index(protein_id).drop_duplicates().sort_values('clumpsptm_pval')

        site_collapse_df = df.reset_index()[
            [protein_id,site_id,'acc_res','acc_res_i',
             'uniprot_res_i','uniprot_res','pdb_res_i','pdb_res','value']
            ].groupby(protein_id).agg(list)

        acc_df = acc_df.join(site_collapse_df)

        # Post-processing
        acc_df['clumpsptm_input_n'] = acc_df['clumpsptm_input'].apply(lambda x: len(x.split("+")))
        acc_df['clumpsptm_sample_n'] = acc_df['clumpsptm_sample'].apply(lambda x: len(x.split("+")))

        # Erroneous
        acc_df = acc_df[acc_df['clumpsptm_input_n']!=acc_df['clumpsptm_sample_n']]

        # Protein - Gene
        acc_df["y"] = acc_df.index + " | " + acc_df['geneSymbol']

        # Filter out for min sites
        acc_df = acc_df[acc_df['clumpsptm_input'].apply(lambda x: len(x.split("+")) >= min_sites)]

        # Threshold lower end of p-values
        if thresh_min_pval:
            acc_df.loc[acc_df['clumpsptm_pval']==0,'clumpsptm_pval'] = 1/acc_df['clumpsptm_niter'][acc_df['clumpsptm_pval']==0]/10

        # FDR
        if acc_df.shape[0] == 0:
            return None

        _,acc_df['clumpsptm_fdr'],_,_ = multipletests(acc_df['clumpsptm_pval'].astype(float), method='fdr_bh', alpha=0.1)

        # Corrected FDR
        acc_df = add_corrected_fdr(acc_df)

        return acc_df

    ptm_files = glob.glob(os.path.join(output_dir, "*_ptm_clumpsptm*"))
    phosph_files = glob.glob(os.path.join(output_dir, "*_phosphoproteome_clumpsptm*"))
    acetyl_files = glob.glob(os.path.join(output_dir, "*_acetylome_clumpsptm*"))

    res_df = list()
    if len(ptm_files) > 0:
        clumps_ptm_df = pd.concat([pd.read_parquet(x) for x in ptm_files])
        clumps_ptm_df['clumpsptm_sampler'] = "ptm"
        clumps_ptm_df = _collapse_sites(clumps_ptm_df)
        res_df.append(clumps_ptm_df)

    if len(phosph_files) > 0:
        clumps_phosph_df = pd.concat([pd.read_parquet(x) for x in phosph_files])
        clumps_phosph_df['clumpsptm_sampler'] = "phosphoproteome"
        clumps_phosph_df = _collapse_sites(clumps_phosph_df)
        res_df.append(clumps_phosph_df)

    if len(acetyl_files) > 0:
        clumps_acetyl_df = pd.concat([pd.read_parquet(x) for x in acetyl_files])
        clumps_acetyl_df['clumpsptm_sampler'] = "acetylome"
        clumps_acetyl_df = _collapse_sites(clumps_acetyl_df)
        res_df.append(clumps_acetyl_df)

    if len(res_df) == 0:
        raise ValueError('NO RESULTS FILES FOUND.')

    return pd.concat(res_df)
