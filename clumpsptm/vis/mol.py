# -- import packages: --------------------------------------------------------------------
import os
from typing import Union
import pandas as pd
from tqdm import tqdm

def buildPymol(
    pdb: str,
    chain: str,
    session_name: str,
    phosph_residues: Union[None, list] = None,
    acetyl_residues: Union[None, list] = None,
    chain_color: str = 'palecyan',
    phosph_color: str = 'deeppurple',
    acetyl_color: str = 'hotpink'
    ):
    """
    Create a PyMol session using python API.

    Parameters:
    -----------
    pdb
        pdb protein name
        type: str
    
    chain
        chain of pdb structure
        type: str

    session_name
        file name to save as (will append .pse)
        type: str

    phosph_residues [optional]
        protein residues to highlight in pymol structure
        type: [None, list]

    acetyl_residues [optional]
        protein residues to highlight in pymol structure
        type: [None, list]

    chain_color [optional]
        color of entire protein structure
        type: str
        default: palecyan

    phosph_color [optional]
        color of phosph highlighted spheres
        type: str
        default: deeppurple

    acetyl_color [optional]
        color of acetyl highlighted spheres
        type: str
        default: hotpink

    Returns:
    --------
    None

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
    pymol.cmd.color(chain_color, obj_name)
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

def buildPymol_from_result(entry: pd.Series, out_dir: Union[None, str] = None, include_idx_in_name: bool = False):
    """
    Build Pymol from entry.
    """
    import ast

    # Name of Pymol Session
    session_name = list()
    if include_idx_in_name:
        session_name.append(entry.name)

    session_name.append(entry['geneSymbol'])
    session_name.append(entry['clumpsptm_sampler'])

    if 'id' in entry:
        session_name.append(entry['id'])

    session_name = '_'.join(session_name)

    # Make directory
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        session_name = os.path.join(out_dir, session_name)

    phosph_residues = list()
    acetyl_residues = list()

    if entry['clumpsptm_sampler'] == 'phosphoproteome':
        phosph_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'acetylome':
        acetyl_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'ptm':
        for (res,num) in zip(ast.literal_eval(entry['pdb_res']), entry['clumpsptm_input'].split('+')):
            if res=='K':
                acetyl_residues.append(str(num))
            else:
                phosph_residues.append(str(num))

    if len(acetyl_residues)==0:
        acetyl_residues = None
    if len(phosph_residues)==0:
        phosph_residues = None

    buildPymol(
        entry['pdb'],
        entry['chain'],
        session_name,
        phosph_residues = phosph_residues,
        acetyl_residues = acetyl_residues
    )

def create_pymols_from_result(results_df: pd.DataFrame, out_dir: str = '.', pval_thresh: float = 0.1, include_idx_in_name: bool = False, verbose: bool = False):
    """Create pymol sessions for all results."""
    _df = results_df[results_df['clumpsptm_pval']<pval_thresh]

    if out_dir!='.':
        if verbose:
            print("   * writing out session files to {}".format(out_dir))

    for idx,row in tqdm(_df.iterrows(), total=_df.shape[0]):
        buildPymol_from_result(row, out_dir=out_dir, include_idx_in_name=include_idx_in_name)
