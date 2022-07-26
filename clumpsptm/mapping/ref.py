from typing import Union
import glob
import os
from sys import stdout
import subprocess
from tqdm import tqdm
from collections import defaultdict
from agutil.parallel import parallelize2
import pandas as pd

from .blast import Blast
from ..utils import AMINO_ACID_MAP

def dl_ref(
    out_dir: str,
    files: Union[tuple, None] = None,
    verbose: bool = True
    ):
    """
    Download reference files
    --------------------------
    Downloads references files to provided directory.
       * PDB Seq Res (for Database)
       * RefSeq Latest Identifiers (fasta)

    """
    os.makedirs(out_dir, exist_ok=True)

    if files is None:
        files = (
            "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt",
            "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_protein.faa.gz",
        )

    for file in files:
        if verbose:
            stdout.write("   * Downloading {}...\n".format(file.split('/')[-1]))
        cmd = "wget -P {} {}".format(out_dir, file)
        subprocess.call(cmd, executable='/bin/bash', shell=True)

        if file.endswith('.gz'):
            cmd = "gunzip -c {} > {}".format(os.path.join(out_dir, file.split('/')[-1]),os.path.join(out_dir, file.split('/')[-1].replace('.gz','')))
            subprocess.call(cmd, executable='/bin/bash', shell=True)

def split_fastas(
    fasta: str,
    out_dir: str,
    naming: str = "cptac"
    ):
    """
    Split fastas into individual files.
    """
    os.makedirs(out_dir, exist_ok=True)
    errors = list()

    with open(fasta, 'r') as f:
        data = f.read().split('>')

        for entry in tqdm(data[1:]):
            if naming=="cptac":
                f_name = entry.split(' ')[0]+'.seq'
            elif naming=="uniprot":
                try:
                    f_name = entry.split('|')[1]+'.seq'
                except:
                    errors.append(entry.split("\n")[0])
                    continue

            with open(os.path.join(out_dir, f_name), 'w') as f:
                f.write('>' + entry)

    if naming == "uniprot":
        print("   * {} sequences not parsed".format(len(errors)))
        return errors

def get_blast_hits_with_sifts(
    blasted_files,
    sifts_df,
    filter_only_human = False
    ):
    """
    Takes a list of blasted files, selects top blasted hit cross-references input
    Sifts data-base file.

    Parameters:
    -----------
    blasted_files [ required ]
        type: list

    sifts_db [ required ]
        type: pd.DataFrame

    filter_only_human [ optional ]
        default: False
        type: bool

    Returns:
    --------
    pd.DataFrame
    """
    import pandas as pd

    uniprots = set(sifts_df['SP_PRIMARY'])

    res = list()

    for bf in tqdm(blasted_files):
        bl = Blast(bf)

        try:
            bl.hits.loc[:,'uniprot'] = bl.hits['Hit_def'].apply(lambda x: x.split("|")[1])
            bl.hits.loc[:,'species'] = bl.hits['Hit_def'].apply(lambda x: x.split("|")[2].split(' ')[0].split("_")[1])
            bl.hits.loc[:,'query'] = bl.query

            if filter_only_human:
                bl.hits = bl.hits[bl.hits['species']=='HUMAN']

            hits_in_sifts = bl.hits[bl.hits['uniprot'].isin(uniprots)]

            if hits_in_sifts.shape[0] == 0:
                _res = bl.hits.iloc[0].copy()
                _res['sifts'] = False
            else:
                _res = hits_in_sifts.iloc[0].copy()
                _res['sifts'] = True

            _res['blast'] = True
        except:
            _res = pd.Series()
            _res['query'] = bl.query
            _res['blast'] = False

        res.append(_res)

    return pd.concat(res,1).T.set_index("query")

def get_pdb_match(sifts_s, acc_df, pdbstore):
    """
    For a given sifts entry, this loads and concatenates PDB information.

    Parameters:
    -----------
    sifts_s [ required ]
        type: pd.Series

    acc_df [ required ]
        type: pd.DataFrame

    pdbstore [ required ]
        type: clumpsptm.pdbstore.PdbStore

    Returns:
    --------
    pd.DataFrame
        accession numbers as indices and corresponding PDB res_num & res
    """
    import numpy as np

    pdb_resnames = pdbstore.load_rd(sifts_s["PDB"], sifts_s["CHAIN"])

    def _try_map(x):
        try:
            return AMINO_ACID_MAP[pdb_resnames[x]]
        except:
            return None

    # PDB IDs
    pdb_i = np.unique(list(pdb_resnames.keys()))

    # Uniprot 2 PDB
    md = defaultdict(int,{x+sifts_s['db_first_to']-sifts_s['db_first_from']:x for x in pdb_i})

    _pdb_res = '{}_{}_res'.format(sifts_s['PDB'], sifts_s['CHAIN'])
    _pdb_res_i = '{}_{}_res_i'.format(sifts_s['PDB'], sifts_s['CHAIN'])
    _pdb_res_match = '{}_{}_res_match'.format(sifts_s['PDB'], sifts_s['CHAIN'])

    acc_df[_pdb_res_i] = acc_df['uniprot_res_i'].apply(lambda x: md[x])
    acc_df[_pdb_res] = acc_df[_pdb_res_i].apply(_try_map, 1)
    acc_df[_pdb_res_match] = acc_df['acc_res']==acc_df[_pdb_res]

    return acc_df[[_pdb_res_i,_pdb_res,_pdb_res_match]]

def get_acc_pdb_matches(
    accession_numbers,
    ptm_df,
    sifts_df,
    pdbstore,
    out_dir,
    n_threads
    ):
    """
    For a list of accession_numbers, use multi-threading to curate mapped PDB
    structures using SIFTS database.

    Parameters:
    -----------
    accession_numbers [ required ]
        type: iterable

    ptm_df [ required ]
        type: pd.DataFrame

    sifts_df [ required ]
        type: pd.DataFrame

    pdbstore [ required ]
        type: clumpsptm.pdbstore.PdbStore

    out_dir [ required ]
        type: str

    n_threads [ required ]
        type: int

    Returns:
    --------
    set
        Accession numbers with errors in mapping.
    """
    # Make directories
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "acc_mapped_sites"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "acc_mapped_sifts"), exist_ok=True)

    # Multi-threading function
    @parallelize2(maximum=n_threads)
    def get_acc_pdb_match(acc):
        """Func for matching"""
        try:
            acc_df = ptm_df[ptm_df['accession_number']==acc].copy().sort_values("uniprot_res_i")[[
                'uniprot',
                'acc_res_i',
                'acc_res_idx',
                'acc_res',
                'uniprot_res_i',
                'uniprot_res',
                'uniprot_match'
            ]]

            uniprot = acc_df.iloc[0]['uniprot']

            sifts_filt_df = sifts_df[sifts_df['SP_PRIMARY']==uniprot].copy().sort_values(['SP_BEG','len'], ascending=[True,False])
            sifts_filt_df = sifts_filt_df.drop_duplicates(subset=['PDB','CHAIN'])

            resd = list()

            for idx,row in sifts_filt_df.iterrows():
                try:
                    _resd = get_pdb_match(row, acc_df, pdbstore)
                    resd.append(_resd)
                    sifts_filt_df.loc[idx,'percent_match'] = 100 * sum(_resd.iloc[:,-1]) / _resd.shape[0]
                except:
                    print("Error loading: {}-{}".format(row["PDB"], row["CHAIN"]))

            resd = pd.concat(resd,1)
            sifts_filt_df['accession_number'] = acc
            sifts_filt_df = sifts_filt_df.astype(str)

            resd.to_parquet(os.path.join(out_dir, "acc_mapped_sites", "{}.parquet".format(acc)))
            sifts_filt_df.to_parquet(os.path.join(out_dir, "acc_mapped_sifts", "{}.parquet".format(acc)))

            return None
        except:
            return acc

    tmp = [get_acc_pdb_match(acc) for acc in accession_numbers]
    err = {callback() for callback in tmp}

    return err
