# CLUMPS-PTM

An algorithm for identifying 3D clusters ("clumps") of post-translational modifications (PTMs). Developed for the Clinical Proteomic Tumor Atlas Consortium ([CPTAC](https://proteomics.cancer.gov/programs/cptac)). Full project repoistory for pan-cancer project can be found [here](https://github.com/getzlab/CPTAC_PanCan_2021).

__Author__: Shankara Anand

__Email__: sanand@broadinstitute.org

_Requires Python 3.6.0 or higher._

## Installation

##### PIP

`pip3 install clumpsptm`

or

##### Git Clone

* `git clone https://github.com/broadinstitute/getzlab-SignatureAnalyzer.git`
* `cd CLUMPS-PTM`
* `pip3 install -e .`

## Use

CLUMPS-PTM has 3 general phases of analysis:
1. __Mapping__: taking input PTM proteomic data and mapping them onto PDB structural data
2. __CLUMPS__: running the algorithm for identifying statistically significant clustering of PTM sites
3. __Post-Processing__: post-processing (FDR correction) \& visualization in pymol
