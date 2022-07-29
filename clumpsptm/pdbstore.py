import importlib,sys
import os
import pandas as pd
import numpy as np
from collections import defaultdict
from typing import Union
from scipy.spatial.distance import euclidean
from .utils import AMINO_ACID_MAP, gunzipper
from Bio import pairwise2
import subprocess
import glob
from agutil.parallel import parallelize2

class PdbStore(object):
    """
    Protein Data Bank Direcotry Store
    ------------------
    This object is used to keep track of protein data bank files
    that are downloaded in a reference directory.
    """
    def __init__(self, path):
        """
        Args:
            * path: PDB directory

        """
        self.pdb_root_dir = path
        self.pdb_dir = os.path.join(path, 'ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/')

    def __str__(self):
        return "PDBStore\n   * {} PDB files downloaded\n   * Directory Path: {}".format(len(self.downloaded_pdbs), self.pdb_root_dir)

    @property
    def downloaded_pdbs(self):
        """
        Get Downloaded Structures.
        ------------------
        Iterates through pdb structure direcotry to find what
        pdbs are already downloaded. Returns a set of downloaded PDBs.
        """
        try:
            return {
                    i[3:7]
                    for j in os.listdir(self.pdb_dir)
                    for i in os.listdir(os.path.join(self.pdb_dir, j))
                    if i.endswith('.ent.gz') and os.path.getsize(os.path.join(self.pdb_dir, j, i)) > 0
                   }
        except:
            return set()

    def missing_pdbs(self, pdbs: set):
        """
        Missing PDBs
        ------------------
        Set difference of input pdbs and all downloaded pdbs.

        Args:
            * pdbs: set of PDBs
        """
        return pdbs - self.downloaded_pdbs

    def download_missing_pdbs(self, pdbs_to_download: Union[set,list], n_threads=15):
        """
        Download Missing PDBs
        ------------------
        Args:
            * pdbs_to_download: a list or set of PDBs to download
            * n_threads: number of threads to use for downloads
        """
        if isinstance(pdbs_to_download, list):
            pdbs_to_download = set(pdbs_to_download)

        out_dir = self.pdb_root_dir

        @parallelize2(maximum=n_threads)
        def dl_pdb(pdb):
            """Download pdb."""
            try:
                cmd = "wget -qr ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/{}/pdb{}.ent.gz -P {}".format(pdb[1:3], pdb, out_dir)
                output = subprocess.check_call(cmd, executable='/bin/bash', shell=True)
                return None
            except:
                return pdb

        print("   * Downloading {} pdbs using {} threads".format(len(pdbs_to_download), n_threads))

        tmp = [dl_pdb(pdb) for pdb in pdbs_to_download]
        pdb_err = {callback() for callback in tmp}

        print("   * Downloaded {} / {} successfully.".format(len(pdbs_to_download-pdb_err) , len(pdbs_to_download)))

        return pdb_err

    def load(self, pdb, chain=None):
        """
        Load residues - amino acids.

        Returns:
            * PDBStream
        """
        from prody import parsePDBStream

        pdb_file = os.path.join(self.pdb_dir, pdb[1:3], "pdb{}.ent.gz".format(pdb))

        with gunzipper(pdb_file) as pfile:
            aa = parsePDBStream(pfile, chain=chain)

        return aa

    def load_header(self, pdb, chain=None):
        """
        Load PDB header

        Returns:
            * dict
        """
        from prody import parsePDBHeader

        pdb_file = os.path.join(self.pdb_dir, pdb[1:3], "pdb{}.ent.gz".format(pdb))
        return parsePDBHeader(pdb_file)

    def load_rd(self, pdb, chain=None):
        """
        Load residue dict
        """
        aa = self.load(pdb,chain)
        res_map = dict(zip(aa.getResnums(),aa.getResnames()))
        return {k:v for k,v in res_map.items() if v in AMINO_ACID_MAP}

    def load_dm(self, pdb, chain=None, point='centroid', return_centroid_coord=False):
        """
        Load distance matrix for pdb.
        """
        from scipy.spatial.distance import euclidean
        import numpy as np

        aa = self.load(pdb,chain=chain)

        xx = aa.getResnums()
        yy = aa.getCoords()
        zz = aa.getResnames()

        pdb_resids = {}
        for i in range(len(xx)):
            if zz[i] in AMINO_ACID_MAP:
                pdb_resids[xx[i]] = True

        mapped_pdb_to_aa = defaultdict(set)
        for idx,resnum in enumerate(xx):
            if resnum in pdb_resids:
                mapped_pdb_to_aa[resnum].add(zz[idx])

        coords = {}
        for i in range(len(xx)):
            if xx[i] not in pdb_resids:
                continue
            if xx[i] not in coords:
                coords[xx[i]] = []
            coords[xx[i]].append(yy[i])  ## add coordinates of an atom belonging to this residue

        ## Euclidean distance matrix
        D = []
        for i in range(len(pdb_resids)):
            D.append(np.zeros(i, dtype=np.float32))

        centroids = {}
        for k in coords:
            centroids[k] = np.mean(np.array(coords[k]), 0)

        co = [centroids[i] for i in pdb_resids]  ## pdb residue coordinates

        for i in range(len(pdb_resids)):
            for j in range(i):
                D[i][j] = euclidean(co[i], co[j])

        if return_centroid_coord:
            return (D, pdb_resids, mapped_pdb_to_aa, co)
        else:
            return (D, pdb_resids, mapped_pdb_to_aa)

class AlphaStore(object):
    """
    Alpha Fold Data Bank Directory Store
    ------------------
    This object is used to keep track of protein data bank files
    that are downloaded in a reference directory.
    """
    def __init__(self, path):
        """
        Args:
            * path: PDB directory

        """
        self.pdb_dir = path
        self.uniprot_dict = {x.split("-")[1]:x for x in glob.glob(os.path.join(self.pdb_dir, "*.pdb.gz"))}

        try:
            self.tar_file = glob.glob(os.path.join(self.pdb_dir, "*.tar"))[0]
            self.version = self.tar_file.split(".tar")[0].split("_")[-1]
        except:
            pass

    def __str__(self):
        return "AlphaStore\n   * {} PDB files downloaded\n   * Directory Path: {}".format(len(self.uniprots), self.pdb_dir)

    @property
    def uniprots(self):
        """
        Return downloaded uniprot models.
        """
        return list(self.uniprot_dict.keys())

    def load(self, uniprot):
        """
        Load residues - amino acids.

        Returns:
            * PDBStream
        """
        from prody import parsePDBStream

        pdb_file = self.uniprot_dict[uniprot]

        with gunzipper(pdb_file) as pfile:
            aa = parsePDBStream(pfile)

        return aa

    def load_header(self, uniprot):
        """
        Load PDB header

        Returns:
            * dict
        """
        from prody import parsePDBHeader

        pdb_file = self.uniprot_dict[uniprot]
        return parsePDBHeader(pdb_file)

    def scrape_header(self, uniprot):
        """
        Scrape PDB header. Just gets alignment information for pdb.
        """
        import re
        import gzip

        header = dict()
        with gzip.open(self.uniprot_dict[uniprot], mode='rt') as pfile:
            for line in pfile.readlines():
                if line[:5] == "DBREF":
                    l = line.strip()
                    l = re.split("\s+", l, flags=re.UNICODE)
                    header[uniprot]["db_first_from"] = int(l[3])
                    header[uniprot]["db_first_to"] = int(l[8])
                    header[uniprot]["db_accession"] = l[6]
                    header[uniprot]["db_ref"] = l[5].replace("UNP","UniProt")

        return header

    def load_rd(self, uniprot):
        """
        Load residue dict
        """
        aa = self.load(uniprot)
        res_map = dict(zip(aa.getResnums(),aa.getResnames()))
        return {k:v for k,v in res_map.items() if v in AMINO_ACID_MAP}

    def load_dm(self, uniprot, point='centroid', return_centroid_coord=False):
        """
        Load distance matrix for pdb.
        """
        from scipy.spatial.distance import euclidean
        import numpy as np

        aa = self.load(uniprot)

        xx = aa.getResnums()
        yy = aa.getCoords()
        zz = aa.getResnames()

        # Model confidence
        bb = aa.getBetas()

        pdb_resids = {}
        for i in range(len(xx)):
            if zz[i] in AMINO_ACID_MAP:
                pdb_resids[xx[i]] = True

        mapped_pdb_to_aa = defaultdict(set)
        model_confidence = defaultdict(float)

        for idx,resnum in enumerate(xx):
            if resnum in pdb_resids:
                mapped_pdb_to_aa[resnum].add(zz[idx])
                model_confidence[resnum] = bb[idx]

        coords = {}
        for i in range(len(xx)):
            if xx[i] not in pdb_resids:
                continue
            if xx[i] not in coords:
                coords[xx[i]] = []
            coords[xx[i]].append(yy[i])  ## add coordinates of an atom belonging to this residue

        ## Euclidean distance matrix
        D = []
        for i in range(len(pdb_resids)):
            D.append(np.zeros(i, dtype=np.float32))

        centroids = {}
        for k in coords:
            centroids[k] = np.mean(np.array(coords[k]), 0)

        co = [centroids[i] for i in pdb_resids]  ## pdb residue coordinates

        for i in range(len(pdb_resids)):
            for j in range(i):
                D[i][j] = euclidean(co[i], co[j])

        if return_centroid_coord:
            return (D, pdb_resids, mapped_pdb_to_aa, model_confidence, co)
        else:
            return (D, pdb_resids, mapped_pdb_to_aa, model_confidence)

class AccessionNo(object):
    """
    Accession Number Object.
    ---------------
    Stores Accession Number -> PDB mapping.
    """
    def __init__(
        self,
        accession_number,
        blast_dir="ref/refseq_blasted",
        pdb_dir="../../../getzlab-CLUMPS2/clumps/db/ref/pdbs"
    ):
        """
        Args:
            * x
        """
        xpo_param=(3,4.5,6,8,10)

        self.accession_number = accession_number
        self.blast_dir = blast_dir
        self.pdb_dir = pdb_dir

        # Get Blasted PDB
        self.blast = clumpsptm.mp.Blast(os.path.join(self.blast_dir,"{}.blasted.seq.gz".format(self.accession_number)))
        self.hit = self.blast.hits.iloc[0]
        self.pdb,self.chain = self.hit['Hit_def'].split(' ')[0].split("_")

        # Get PDB Info
        pdbstore = clumpsptm.PdbStore(self.pdb_dir)
        self.D,self.x,self.pdb_resnames = pdbstore.load_dm(self.pdb,self.chain)
        self.DDt = clumpsptm.transform_dx(self.D, xpo_param)

        # Align PDB Fasta to PDB SEQATOMS
        self.pdb_structure_idx = np.arange(max(list(self.pdb_resnames.keys())))+1
        self.pdb_structure_res = np.array(["."]*pdb_idx.shape[0])

        for idx in self.pdb_structure_idx-1:
            pdb_i = self.pdb_structure_idx[idx]
            try:
                self.pdb_structure_res[idx] = AMINO_ACID_MAP[list(self.pdb_resnames[pdb_i])[0]]
            except:
                pass

        self.pdb_structure_res = ''.join(self.pdb_structure_res)
        self.pdb_fasta_alignment = pairwise2.align.globalms(
            self.hit['Hsp_hseq'],
            self.pdb_structure_res.replace(".",""),1, -1, -.5, -.1
        )[0]
        self.offset = self.get_init_pdb_struct_i(self.pdb_structure_res) - self.get_alignment_offset(self.pdb_fasta_alignment[1])

    def get_init_pdb_struct_i(self, pdb_structure_res):
        """
        Get index of first registered atom.
        """
        for idx,a in enumerate(pdb_structure_res):
            if a!='.':
                return idx+1

    def get_alignment_offset(self, align):
        """
        Get offset from Fasta via alignment
        """
        c = 0
        for a in align:
            if a=='-':
                c+=1
            else:
                return c
