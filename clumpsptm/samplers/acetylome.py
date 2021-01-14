import numpy as np

class AcetylSampler(object):
    """
    Samples only Lysine residues.
    --------------------------
    Takes in pdb residue id (pr) and a mapping default dictionary.
    Samples only from Lysines.
    """
    def __init__(self, pdb_resid_to_d_map, pdb_resnames, res='LYS'):
        self.resid_to_d_map = pdb_resid_to_d_map
        self.mapping = pdb_resnames

        self.lysines = list()

        for idx in pdb_resnames.keys():
            try:
                if list(pdb_resnames[idx])[0]==res:
                    self.lysines.append(idx)
            except:
                pass

    def sample(self, mi):
        sorted_indices = sorted(random.sample(self.lysines, len(mi)))
        random_permuts = np.random.permutation(len(mi))
        return [self.resid_to_d_map[x] for x in sorted_indices], random_permuts
