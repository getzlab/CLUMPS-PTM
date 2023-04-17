# CLUMPS-PTM

An algorithm for identifying 3D clusters ("clumps") of post-translational modifications (PTMs). Developed for the Clinical Proteomic Tumor Atlas Consortium ([CPTAC](https://proteomics.cancer.gov/programs/cptac)). Full project repoistory for pan-cancer project can be found [here](https://github.com/getzlab/CPTAC_PanCan_2021).

*This tool is still in early development, if you are interested in collaborating, please reach out to the email below.*

__Author__: Shankara Anand

__Email__: sanand@broadinstitute.org

_Requires Python 3.6.0 or higher._

## Installation

#### PIP

`pip3 install clumps-ptm`

or

#### Git Clone

```
git clone git@github.com:getzlab/CLUMPS-PTM.git
cd CLUMPS-PTM
pip3 install -e .
```

## Use

CLUMPS-PTM has 3 general phases of analysis:
1. __Mapping__: taking input PTM proteomic data and mapping them onto PDB structural data.

  Mapping relies on the source data and involves programmatic calling of `blastp+` depending on the source data-base to map to UNIPROT and ultimately PDB structures. An example notebook that walks through the mapping and demonstrates use of `clumps-ptm` API for running these steps programmatically can be found [here](https://github.com/getzlab/CLUMPS-PTM/blob/main/examples/CPTAC_Mapping_Workflow.ipynb). Once the mapping is performed once for a new data-set, the mapping file is used as the `--maps` flag in `clumpsptm` command (below). Example workflow using AlphaFold structures may be found [here](https://github.com/getzlab/CPTAC_PanCan_2021/blob/master/clumpsptm_analysis/pancan_alphafold_clumpsptm/01_CPTAC_pdb_workflow.ipynb).

2. __CLUMPS__: running the algorithm for identifying statistically significant clustering of PTM sites.

CLUMPS-PTM was designed for either 1) differential expression proteomic data or 2) sample frequency smoothed weights.

[1] Differential expression: due to the nature of drop-out in Mass-Spectrometry data, we opt for using broad changes in PTM levels across sample groups to interrogate "clumping" of modifications. We use output from Limma-Voom differential expression, and set `n` (input weights) to `logFC x -log10FDR`.

[2] Sample frequency: compute smoothed sample frequency using Sigmoidal Hill function (similiar to the original CLUMPS based on missense mutations). For example:

```
def compute_freq(df, theta=12, m=3):
    """Compute frequency file."""
    freq = pd.DataFrame(df.notna().sum(1).sort_values(ascending=False), columns=['raw'])
    freq = freq[freq['raw']>0]
    freq['n'] = np.power(freq['raw'], m) / (np.power(theta, m) + np.power(freq['raw'],m))
    
    return freq
 ```
  
Example input:
  
|                                    |        n | feature         | id   |
|:-----------------------------------|---------:|:----------------|:-----|
| NP_001128690.1_S51s_1_1_51_51      | 0.85 | phosphoproteome | GroupA |
| NP_001309339.1_S211s_1_1_211_211   | 0.22 | phosphoproteome | GroupA |
| NP_056988.3_S214s_1_1_214_214      | 0.1 | phosphoproteome | GroupA |
| NP_056988.3_S164s_1_1_164_164      | 0.0.43 | phosphoproteome | GroupA |

*Notes:*
* `feature` is required to indicate the type of PTM sampler to use
* `id` is required to indicate the subset of group to use
* if `--subset` is set to positive or negative, the `n` column will be filtered to run clumps for either up-regulated or down-regulated sites when using differential expression

#### Command Line

```{python}
usage: clumpsptm [-h] -i INPUT -m MAPS -w WEIGHT -s PDBSTORE [-o OUTPUT_DIR]
                 [-x XPO] [--threads THREADS] [-v]
                 [-f [FEATURES [FEATURES ...]]] [-g GROUPING] [-q]
                 [--min_sites MIN_SITES] [--subset {positive,negative}]
                 [--protein_id PROTEIN_ID] [--site_id SITE_ID] [--alphafold]
                 [--alphafold_threshold ALPHAFOLD_THRESHOLD]

Run CLUMPS-PTM.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        <Required> Input file.
  -m MAPS, --maps MAPS  <Required> Mapping with index as indices that overlap
                        input.
  -w WEIGHT, --weight WEIGHT
                        <Required> Weighting for CLUMPS-PTM (ex. logFC).
  -s PDBSTORE, --pdbstore PDBSTORE
                        <Required> path to PDBStore directory.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory.
  -x XPO, --xpo XPO     Soft threshold parameter for truncated Gaussian.
  --threads THREADS     Number of threads for sampling.
  -v, --verbose         Verbosity.
  -f [FEATURES [FEATURES ...]], --features [FEATURES [FEATURES ...]]
                        Assays to subset for.
  -g GROUPING, --grouping GROUPING
                        DE group to use.
  -q, --use_only_significant_sites
                        Only use significant sites for CLUMPS-PTM.
  --min_sites MIN_SITES
                        Minimum number of sites.
  --subset {positive,negative}
                        Subset sites.
  --protein_id PROTEIN_ID
                        Unique protein id in input.
  --site_id SITE_ID     Unique site id in input.
  --alphafold           Run using alphafold structures.
  --alphafold_threshold ALPHAFOLD_THRESHOLD
                        Threshold confidence level for alphafold sites.
                        
```

Example scripts used for the CPTAC project where this is run using differential expression results may be found [here](https://github.com/getzlab/CPTAC_PanCan_PTM_2023/blob/master/CLUMPS-PTM/run.sh).

FDR correction is [done](https://github.com/getzlab/CLUMPS-PTM/blob/5713948a0398372ffb1bcadff122b93eedc90b76/clumpsptm/utils.py#L40) when running the tool from command line.

3. __Visualization__: both Pymol structures + p-value hits.

```
import clumpsptm

# Create p-value dotplot
clumpsptm.vis.dotplot(...)

# Create pymol session with CLUMPS-PTM hit
clumpsptm.vis.buildPymol(...)
```
