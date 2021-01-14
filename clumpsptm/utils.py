import os
from gzip import GzipFile
from prody import confProDy
import contextlib
import tempfile
import subprocess

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
