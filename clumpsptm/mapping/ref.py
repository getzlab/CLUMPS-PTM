from typing import Union
import glob
import os
from sys import stdout
import subprocess
from tqdm import tqdm

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
    out_dir: str
    ):
    """
    Split fastas into individual files.
    """
    os.makedirs(out_dir, exist_ok=True)

    with open(fasta, 'r') as f:
        data = f.read().split('>')

        for entry in tqdm(data[1:]):
            with open(os.path.join(out_dir, entry.split(' ')[0]+'.seq'), 'w') as f:
                f.write('>' + entry)
