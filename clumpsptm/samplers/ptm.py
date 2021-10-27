import numpy as np
import random

class PTMSampler(object):
    """
    Samples only modifiable residues.
    --------------------------
    Takes in pdb residue id (pr) and a mapping default dictionary.
    Samples only from residues that can be acetylated or phosphorylated.
    """
    def __init__(self, pdb_resid_to_d_map, pdb_resnames, res=['LYS','SER','THR','TYR']):
        self.resid_to_d_map = pdb_resid_to_d_map
        self.mapping = pdb_resnames
        self.ptm = list()

        for idx in pdb_resnames.keys():
            try:
                if list(pdb_resnames[idx])[0] in res:
                    self.ptm.append(idx)
            except:
                pass

    def sample(self, mi):
        sorted_indices = sorted(random.sample(self.ptm, len(mi)))
        random_permuts = np.random.permutation(len(mi))
        return [self.resid_to_d_map[x] for x in sorted_indices], random_permuts
