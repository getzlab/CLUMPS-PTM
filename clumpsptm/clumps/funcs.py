import sys
import os
import pandas as pd
import numpy as np

def transform_dx(D, xpo_param):
    """
    Transform distance matrix.
    --------------------------

    """
    import numpy as np

    DDt = []
    for soft_thresh_idx in xpo_param:
        den = 2.0 * soft_thresh_idx**2
        m = []
        for i in range(len(D)):
            mrow = np.zeros(i, dtype=np.float32)
            for j in range(i):
                mrow[j] = np.exp(-(D[i][j]**2)/den)
            m.append(mrow)
        DDt.append(m)

    return DDt

def wap(mut_indices, mvcorr, Mmv, DDt):
    """
    WAP Score
    --------------------------
    Compute WAP score to summarize pairwise distances between
    mutated residues in 3D protein structure.
    """
    import numpy as np

    s = np.zeros(len(DDt), np.float64)
    for mat in range(len(DDt)):
        d = DDt[mat]
        for i in range(len(mut_indices)):
            dcol = d[mut_indices[i]]
            for j in range(i):
                s[mat] += Mmv[mvcorr[i]][mvcorr[j]] * dcol[mut_indices[j]]
    return s

def init_alg(clumps_input_df, DDt):
    """
    Initialize CLUMPS algorithm.
    """
    mi = np.array(clumps_input_df['d'])
    mv = np.array(clumps_input_df['value'])

    # Compute matrix
    Mmv = []
    mvcorr = range(len(mv))

    for i in range(len(mi)):
        mrow = np.zeros(len(mi), np.float64)
        for j in range(len(mi)):
            mrow[j] = mv[i]*mv[j]
        Mmv.append(mrow)

    # Compute WAP score
    wap_obs = wap(mi, mvcorr, Mmv, DDt)

    # Init Params
    init_dict = {}
    init_dict['mi'] = mi
    init_dict['mv'] = mv
    init_dict['Mmv'] = Mmv
    init_dict['mvcorr'] = mvcorr
    init_dict['wap_obs'] = wap_obs

    return init_dict
