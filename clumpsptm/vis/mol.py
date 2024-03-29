# -- import packages: --------------------------------------------------------------------
import os
from typing import Union
import pandas as pd
from tqdm import tqdm
import numpy as np

def buildPymol(
    pdb: str,
    chain: str,
    session_name: str,
    phosph_residues: Union[None, list] = None,
    acetyl_residues: Union[None, list] = None,
    ubiq_residues: Union[None, list] = None,
    chain_color: str = 'palecyan',
    phosph_color: str = 'deeppurple',
    acetyl_color: str = 'hotpink',
    ubiq_color: str = 'red'
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

    if ubiq_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(ubiq_residues)))
        pymol.cmd.color(ubiq_color, 'resi {}'.format('+'.join(ubiq_residues)))

    pymol.cmd.save(session_name+'.pse')
    os.remove(temp_name)

def buildAlphaPymol(
    uniprot: str,
    alphaStore,
    session_name: str,
    phosph_residues: Union[None, list] = None,
    acetyl_residues: Union[None, list] = None,
    ubiq_residues: Union[None, list] = None,
    chain_color: str = 'palecyan',
    phosph_color: str = 'deeppurple',
    acetyl_color: str = 'hotpink',
    ubiq_color: str = 'red',
    color_chain_by_confidence: bool = False,
    ):
    """
    Create a PyMol session using python API.

    Parameters:
    -----------
    uniprot
        uniprot ID
        type: str

    alphaStore
        alphaStore with downloaded structures
        type: AlphaStore

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

    pymol.finish_launching()
    pymol.cmd.delete("all")
    pymol.cmd.load(alphaStore.uniprot_dict[uniprot])
    obj_name = alphaStore.uniprot_dict[uniprot].split("/")[-1].split(".")[0]

    # Stylings
    pymol.cmd.set('ray_trace_fog','0')
    pymol.cmd.set('ray_shadows','0')
    pymol.cmd.set('depth_cue','0')
    pymol.cmd.bg_color('white')
    pymol.cmd.set('antialias','4')
    pymol.cmd.set('cartoon_transparency','0')
    pymol.cmd.set('transparency','0.7')
    pymol.cmd.set('ray_trace_mode','3')

    # Color chains by confidnece based on Alphafold pLDDT
    if color_chain_by_confidence:
        _,_,_,model_confidence = alphaStore.load_dm(uniprot)

        def conf_lvl(x):
            """
            Confidence Level
                Very high (pLDDT > 90)
                High (90 > pLDDT > 70)
                Low (70 > pLDDT > 50)
                Very low (pLDDT < 50)
            """
            if x>= 90:
                return "very high"
            elif x>= 70:
                return "high"
            elif x>= 50:
                return "low"
            else:
                return "very low"

        # Colors
        colors = {"very high":"tv_blue","high":"lightblue","low":"yelloworange","very low":"orange"}

        _df = pd.DataFrame((model_confidence.keys(),model_confidence.values())).T
        _df['conf'] = _df.loc[:,1].apply(conf_lvl)
        _df['colors'] = _df['conf'].apply(lambda x: colors[x])

        for color in np.unique(_df['colors']):
            pymol.cmd.color(color, 'resi {}'.format('+'.join(_df[_df['colors']==color][0].astype(int).astype(str))))

    else:
        pymol.cmd.color(chain_color, obj_name)

    pymol.cmd.show('surface', obj_name)
    pymol.cmd.remove('solvent')

    if phosph_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(phosph_residues)))
        pymol.cmd.color(phosph_color, 'resi {}'.format('+'.join(phosph_residues)))

    if acetyl_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(acetyl_residues)))
        pymol.cmd.color(acetyl_color, 'resi {}'.format('+'.join(acetyl_residues)))

    if ubiq_residues is not None:
        pymol.cmd.show('spheres', 'resi {}'.format('+'.join(ubiq_residues)))
        pymol.cmd.color(ubiq_color, 'resi {}'.format('+'.join(ubiq_residues)))

    pymol.cmd.save(session_name+'.pse')

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
    ubiq_residues = list()

    if entry['clumpsptm_sampler'] == 'phosphoproteome':
        phosph_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'acetylome':
        acetyl_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'ubiquitylome':
        ubiq_residues = entry['clumpsptm_input'].split('+')
    elif entry['clumpsptm_sampler'] == 'ptm':
        # Fix for ubiquitylation sites
        for (res,num) in zip(ast.literal_eval(entry['pdb_res']), entry['clumpsptm_input'].split('+')):
            if res=='K':
                acetyl_residues.append(str(num))
            else:
                phosph_residues.append(str(num))

    if len(phosph_residues)==0:
        phosph_residues = None
    if len(acetyl_residues)==0:
        acetyl_residues = None
    if len(ubiq_residues)==0:
        ubiq_residues = None

    buildPymol(
        entry['pdb'],
        entry['chain'],
        session_name,
        phosph_residues = phosph_residues,
        acetyl_residues = acetyl_residues,
        ubiq_residues = ubiq_residues
    )

def create_pymols_from_result(results_df: pd.DataFrame, out_dir: str = '.', pval_thresh: float = 0.1, include_idx_in_name: bool = False, verbose: bool = False):
    """Create pymol sessions for all results."""
    _df = results_df[results_df['clumpsptm_pval']<pval_thresh]

    if out_dir!='.':
        if verbose:
            print("   * writing out session files to {}".format(out_dir))

    for idx,row in tqdm(_df.iterrows(), total=_df.shape[0]):
        buildPymol_from_result(row, out_dir=out_dir, include_idx_in_name=include_idx_in_name)
