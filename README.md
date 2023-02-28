# CLUMPS-PTM

An algorithm for identifying 3D clusters ("clumps") of post-translational modifications (PTMs). Developed for the Clinical Proteomic Tumor Atlas Consortium ([CPTAC](https://proteomics.cancer.gov/programs/cptac)). Full project repoistory for pan-cancer project can be found [here](https://github.com/getzlab/CPTAC_PanCan_2021).

__Author__: Shankara Anand

__Email__: sanand@broadinstitute.org

_Requires Python 3.6.0 or higher._

## Installation

##### PIP

`pip3 install clumps-ptm`

or

##### Git Clone

```
git clone git@github.com:getzlab/CLUMPS-PTM.git
cd CLUMPS-PTM
pip3 install -e .
```

## Use

CLUMPS-PTM has 3 general phases of analysis:
1. __Mapping__: taking input PTM proteomic data and mapping them onto PDB structural data.

  Mapping relies on the source data and involves programmatic calling of `blastp+` depending on the source data-base to map to UNIPROT and ultimately PDB structures. An example notebook that walks through the mapping and demonstrates use of `clumps-ptm` API for running these steps programmatically can be found [here](https://github.com/getzlab/CLUMPS-PTM/blob/main/examples/CPTAC_Mapping_Workflow.ipynb). Once the mapping is performed once for a new data-set, the mapping file is used as the `--maps` flag in `clumpsptm` command (below).

2. __CLUMPS__: running the algorithm for identifying statistically significant clustering of PTM sites.

  CLUMPS-PTM was designed for use with differential expression proteomic data. Due to the nature of drop-out in Mass-Spectrometry data, we opt for using broad changes in PTM levels across sample groups to interrogate "clumping" of modifications. Thus, the input requires out-put from Limma-Voom differential expression.

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

3. __Post-Processing__: post-processing (FDR correction) \& visualization in Pymol.
