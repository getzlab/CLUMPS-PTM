import os
import sys
from typing import Union
import pandas as pd
from tqdm import tqdm

def buildPymol(
    pdb: str,
    chain: str,
    session_name: str,
    phosph_residues: Union[None, list] = None,
    acetyl_residues: Union[None, list] = None,
    chain_color: str = 'palyecyan',
    phosph_color: str = 'deeppurple',
    acetyl_color: str = 'hotpink'
    ):
    """
    Build Pymol Session
    """
    import __main__
    import tempfile
    import os

    # Quiet and no GUI
    __main__.pymol_argv = [ 'pymol', '-qc']

    import pymol

    temp_name = next(tempfile._get_candidate_names())+'.cif'
    obj_name = "{}{}".format(pdb, chain)

    pymol.finish_launching()

    pymol.cmd.delete("all")
    pymol.cmd.fetch(obj_name, file=temp_name)

    # Stylings
    pymol.cmd.set('ray_trace_fog','0')
    pymol.cmd.set('ray_shadows','0')
    pymol.cmd.set('depth_cue','0')
    pymol.cmd.bg_color('white')
    pymol.cmd.set('antialias','4')
    pymol.cmd.set('cartoon_transparency','0')
    pymol.cmd.set('transparency','0.7')
    pymol.cmd.set('ray_trace_mode','3')
    pymol.cmd.color('palecyan', obj_name)
    pymol.cmd.show('surface', obj_name)
    pymol.cmd.remove('solvent')

    if phosph_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(phosph_residues)))
        pymol.cmd.color(phosph_color, 'resi {}'.format('+'.join(phosph_residues)))

    if acetyl_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(acetyl_residues)))
        pymol.cmd.color(acetyl_color, 'resi {}'.format('+'.join(acetyl_residues)))

    pymol.cmd.save(session_name+'.pse')
    os.remove(temp_name)

def buildPymol_from_result(entry: pd.Series, out_dir: Union[None, str] = None):
    """
    Build Pymol from entry.
    """
    if 'id' in entry:
        session_name = entry['geneSymbol']+'_'+entry['clumpsptm_sampler'] + entry['id']
    else:
        session_name = entry['geneSymbol']+'_'+entry['clumpsptm_sampler']

    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        session_name = os.path.join(out_dir, session_name)

    phosph_residues = None
    acetyl_residues = None

    if entry['clumpsptm_sampler'] == 'phosphoproteome':
        phosph_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'acetylome':
        acetyl_residues = entry['clumpsptm_input'].split('+')

    buildPymol(
        entry['pdb'],
        entry['chain'],
        session_name,
        phosph_residues = phosph_residues,
        acetyl_residues = acetyl_residues
    )

def create_pymols_from_result(results_df: pd.DataFrame, out_dir: str = '.', pval_thresh: float = 0.1):
    """Create pymol sessions for all results."""
    _df = results_df[results_df['clumpsptm_pval']<pval_thresh]

    if out_dir is not '.':
        print("   * writing out session files to {}".format(out_dir))

    for idx,row in tqdm(_df.iterrows(), total=_df.shape[0]):
        buildPymol_from_result(row, out_dir=out_dir)
