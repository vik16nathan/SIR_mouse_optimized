# SIR_mouse_optimized

Code for simulating the amount of alpha synuclein accumulation/brain atrophy originating from an "epicentre."  The model optimizes parameters to simulate the infection and spreading of endogenous aSyn (coded by the Snca gene) across the structural connectome into mesoscale brain regions, predicting output maps of alpha synuclein aggregation or brain atrophy. 

Code is built upon previous work from Shady Rahayel (https://github.com/srahayel/SIR_mouse ), but modified to include the non-uniform clearance of aSyn from each brain region (defined by gene expression), to optimize parameters that best fit our output pathology maps, and to distinguish between aSyn accumulation and downstream atrophy (based on the methods from Zheng et al., 2019). 

## Main script: abm_optuna_general.py

## Full Workflow 

### Filter and package clearance genes into individual model inputs

1. **package_data.R**: provided by Yohan Yee, PhD, to parcellate the coronal and sagittal gene expression patterns for all genes from the Allen Mouse Brain Atlas' in-situ hybridization (ISH) data into the regions in the Allen Brain Ontology. Of particular interest are columns denoting the proportion of "valid" voxels with a measurement in each brain region, as well as the gene expression within each brain region normalized against whole-brain expression. \n
-- output: coronal_expression_data.csv, sagittal_expression_data.csv

2. **filter_parcellated_GE.R**

3. **chop_GE_mr_vv.py**

### Run SIR agent-based model to find optimal parameters and simulate pathology 

4. ***abm_optuna_general.py**: see above
5. ***abm_clearance_genes.py**: allows you to simulate the "I" fraction/atrophy given a specific set of parameters for v, spread_rate, k1, k2, and retro. Output is the .pkl file containing the simulated neurodegeneration data for every timestep of the SIR model, but with a single set of "trial" parameters (no optuna). This is useful for simulating neurodegeneration in each region for the "optimal" parameterization determined above, and then correlating with the empirical pathology.

### Miscellaneous other scripts

## Dependencies
python 3.11.5 with all prerequisite packages installed
R with all prerequisite packages installed
qbatch is not necessary but highly recommended to parallelize job submission (especially when repeating abm_optuna_general.py with thousands of clearance genes and hundreds of trials)

## Main script: abm_optuna_general.py

## Outputs




