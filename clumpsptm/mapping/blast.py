from agutil.parallel import parallelize2
import os
import subprocess
from typing import Union

def blast_sequences(
    db_seq_file: str,
    blast_files: Union[list,str],
    output_dir: Union[None, str] = None,
    db_title: str = "ClumpsRefSeq",
    blast_dir_name = "blast_dir",
    seq_db_name = "seq_db",
    verbose: bool = True,
    n_threads: int = 15,
    ):
    """
    Blast Sequences
    --------------------------
    First, install blast:
        sudo apt-get install ncbi-blast+

    Args:
        * db_seq_file: fasta file to use for database with sequences of inputs
            * ex. PDB sequences, Uniprot
        * blast_files: file(s) of sequences to blast against
        * output_dir: directory to write outputs to
        * db_title: blast database title
        * verbose
        * n_threads: n_threads for blasts

    """
    # Make Blast DIrectory
    os.makedirs(os.path.join(output_dir, blast_dir_name), exist_ok=True)

    # Make Blast Database
    seq_db = os.path.join(output_dir, seq_db_name)
    os.makedirs(seq_db, exist_ok=True)

    # Make Blast Database
    cmd = "makeblastdb -in {} -dbtype 'prot' -out {}/{}/{} -title {}".format(db_seq_file, output_dir, seq_db_name, seq_db_name, db_title)
    db_out = subprocess.run(cmd, executable='/bin/bash', shell=True, stdout=subprocess.PIPE).stdout.decode()

    if verbose: print(db_out)

    @parallelize2(maximum=n_threads)
    def blast(seq):
        """Blast Command."""
        accession_no = seq.split('/')[-1].split('.seq')[0]

        cmd = "blastp -query {} -db {} -out {} -outfmt 5 && gzip -f {} ;".format(
            seq,
            os.path.join(output_dir, seq_db_name, seq_db_name),
            os.path.join(os.path.join(output_dir, blast_dir_name), "{}.blasted.seq".format(accession_no)),
            os.path.join(os.path.join(output_dir, blast_dir_name), "{}.blasted.seq".format(accession_no))
        )
        subprocess.call(cmd, executable='/bin/bash', shell=True)

    if verbose: print("   * Running {} blast sequences using {} threads".format(len(blast_files), n_threads))

    tmp = [blast(s) for s in blast_files]
    results = [callback() for callback in tmp]

class Blast(object):
    """
    Blast results
    ------------------
    This object store blast results.
    """
    def __init__(self, file, eval_thresh=0.001, ident_thresh=0):
        """
        Args:
            * file: blasted filename
            * eval_thresh: threshold for entries
            * ident_thresh: threshold of identity

        """
        from lxml import etree
        import gzip
        import pandas as pd

        self.file_name = file
        self.hits = list()
        self.list = list()
        self.idents = list()

        with gzip.open(file, 'r') as f:
            cur = etree.fromstring(f.read())
            qlen = 0  ## query length

            for i in cur:
                if i.tag == 'BlastOutput_iterations':
                    cur = i
                elif i.tag == 'BlastOutput_query-len':
                    qlen = int(i.text)

            cur = cur[-1]  ## last iteration

            for i in cur:
                if i.tag == 'Iteration_hits':
                    cur = i
                    break

            for hit in cur:
                hit_def = None
                hsps = None

                for i in hit:
                    if i.tag == 'Hit_def':
                        hit_def = i.text
                    elif i.tag == 'Hit_hsps':
                        hsps = i

                for hsp in hsps:
                    di = {}
                    di['Hit_def'] = hit_def

                    for i in hsp:
                        di[i.tag] = i.text

                    if float(di['Hsp_evalue']) > eval_thresh:  ### !!! FILTER !!!
                        continue

                    rel_ident = float(di['Hsp_identity'])/float(di['Hsp_align-len'])

                    if rel_ident < ident_thresh:         ### !!! FILTER !!!
                        continue

                    di['relIdentity'] = rel_ident
                    self.hits.append(di)
                    self.list.append(di['Hit_def'])

        self.hits = pd.DataFrame(self.hits)

    @property
    def nhits(self):
        return len(self.hits)

    @property
    def query(self):
        return self.file_name.split('/')[-1].split('.blasted')[0]

    def __str__(self):
        return "Blastp+ for {} | {} hits".format(self.query, self.nhits)

def load_accession_blast(a_num, unique_sites_df):
    """
    a_num
    """
    an_dict = {}
    an_dict['accession_no'] = a_num

    if scraped_accessions[a_num] is not None:
        an_dict['accession_no_matched'] = scraped_accessions[a_num].split('/')[-1].split('.seq')[0]

        uns = unique_sites_df.loc[a_num,'site_position']

        if isinstance(uns,np.int64):
            an_dict['acetyl_sites'] = np.array([uns,])
        else:
            an_dict['acetyl_sites'] = np.array(uns)

        an_dict['n_acetyl_sites'] = an_dict['acetyl_sites'].shape[0]

        matched = load_blasted_refseq(os.path.join(output_dir, "refseq_blasted", "{}.blasted.seq.gz".format(an_dict['accession_no_matched'])))

        if len(matched) == 0:
            return an_dict

        # Add PDB Chain & Name
        pdb_name = matched[0]['Hit_def'].split(' ')[0]
        an_dict['pdb'] = pdb_name.split('_')[0]
        an_dict['chain'] = pdb_name.split('_')[1]
        an_dict['pdb_desc'] = ' '.join(matched[0]['Hit_def'].split(' ')[1:])
        an_dict['pdb_ident'] = matched[0]['relIdentity']

        # Map Residues
        arr = np.array(list(matched[0]['Hsp_qseq']))

        # Get Overlapping Acetyl-Sites
        overlapping_sites = list()
        res_matches = list()

        for ac in an_dict['acetyl_sites']:
            try:
                res_match = arr[ac-int(matched[0]['Hsp_query-from'])]
            except:
                res_match = None

            if res_match=='K':
                overlapping_sites.append(ac)

            res_matches.append(res_match)

        an_dict['pdb_site_matches'] = overlapping_sites
        an_dict['pdb_n_site_matches']  = len(overlapping_sites)
        an_dict['Hsp_query-from'] = matched[0]['Hsp_query-from']
        an_dict['Hsp_query-to'] = matched[0]['Hsp_query-to']
        an_dict['Hsp_hit-from'] = matched[0]['Hsp_hit-from']
        an_dict['Hsp_hit-to'] = matched[0]['Hsp_hit-to']

        try:
            an_dict['percent_k'] = res_matches.count("K") / (len(an_dict['acetyl_sites'])-res_matches.count(None))
        except:
            an_dict['percent_k'] = None

    return an_dict
