import os
from gzip import GzipFile
from prody import confProDy
import contextlib
import tempfile
import subprocess
import glob
import pandas as pd

AMINO_ACID_MAP = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E', \
         'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',\
         'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','PYX':'C',\
         'SEP':'S','TPO':'T','TYS':'Y','MK8':'L','M3L':'K','DPR':'P','DSN':'S',\
         'ALN':'A','DLY':'K','MVA':'V','MLE':'L','DLE':'L','CSO':'C','PTR':'Y',\
         'BMT':'T','DAL':'A','FAK':'K','MSE':'M','SMF':'A','HYP':'P'}

@contextlib.contextmanager
def gunzipper(gz_file):
    """
    Gunzipper
    --------------------------
    Takes in a file name. Returns a valid file object in an unzipped form. Uses
    context manager to create a temporary file and save into that format. This
    was written for programmatic extraction of PDB files using Prody.

    Inputs:
        gz_file: gzipped file

    Outputs:
        temp_file: temporary file
    """
    with tempfile.NamedTemporaryFile('r', suffix=os.path.splitext(gz_file)[1]) as temp_file:
        subprocess.check_call("gzip -dc {} >> {}".format(gz_file, temp_file.name), executable='/bin/bash', shell=True)
        yield temp_file

def generate_clumpsptm_output(output_dir, min_sites=3):
    """
    Generate CLUMPS-PTM Output file.
    -------------------------------
    Args:
        * path to output_dir
    """
    def _collapse_sites(df):
        """asdfasf."""
        from statsmodels.stats.multitest import multipletests

        acc_df = df[
            ['geneSymbol','uniprot','accession_number','pdb','chain','clumpsptm_sampler',
             'clumpsptm_niter','clumpsptm_pval','clumpsptm_input','clumpsptm_sample']
        ].set_index("accession_number").drop_duplicates().sort_values('clumpsptm_pval')


        site_collapse_df = df.reset_index()[
            ['accession_number','variableSites','acc_res','acc_res_i',
             'uniprot_res_i','uniprot_res','pdb_res_i','pdb_res','value']
            ].groupby("accession_number").agg(list)

        acc_df = acc_df.join(site_collapse_df)

        # Filter out for min sites
        acc_df = acc_df[acc_df['clumpsptm_input'].apply(lambda x: len(x.split("+")) >= min_sites)]

        # FDR
        _,acc_df['clumpsptm_fdr'],_,_ = multipletests(acc_df['clumpsptm_pval'].astype(float), method='fdr_bh', alpha=0.1)

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
