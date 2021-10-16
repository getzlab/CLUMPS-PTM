from typing import Union

from .funcs import wap

def clumps(
    init_dict,
    sam,
    DDt,
    xpo_param: Union[list,tuple],
    seed: int = 1,
    max_rand: int = 1e8,
    use_booster: bool = False,
    booster_interval: int = 1000,
    verbose: bool = False,
    ):
    """
    Run Clumps Algorithm
    ------------
    Returns:
        * p-values
        * wap score aggr
        * n-iterations
    """
    from sys import stdout
    import numpy as np

    def booster():
        """
        Implements the booster algorithm by Getz et al. that saves CPU time.
        Returns False if the randomization should be stopped.
        """
        ret = False
        for i in range(len(xpo_param)):
            s = (rnd - p[i] + 1.0) / ((rnd + 3)*(p[i] + 1))
            if s >= 0.05271:  ## =0.9/(2*1.96)]**2:
                ret = True
                break
        return ret

    if seed is not None:
        np.random.seed(seed)

    p = [0]*len(xpo_param)
    wap_rnd = [0]*len(xpo_param)
    rnd = 0

    for rnd in range(int(max_rand)):
        if rnd%booster_interval == 0 and rnd > 0 and use_booster:
            if verbose:
                stdout.write("\r{}/{}".format(rnd,int(max_rand)))
            if not booster():
                break

        # Randomly permute
        mi,mvcorr = sam.sample(init_dict['mi'])

        # Compute test WAP
        r = wap(mi, mvcorr, init_dict['Mmv'], DDt)

        for rr in range(len(xpo_param)):
            wap_rnd[rr] += r[rr]
            if r[rr] >= init_dict['wap_obs'][rr]:
                p[rr] += 1

    return p, wap_rnd, rnd+1
