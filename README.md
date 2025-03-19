# SIR_mouse_optimized

Code for simulating the amount of alpha synuclein accumulation/brain atrophy originating from an "epicentre."  The model optimizes parameters to simulate the infection and spreading of endogenous aSyn (coded by the Snca gene) across the structural connectome into mesoscale brain regions, predicting output maps of alpha synuclein aggregation or brain atrophy. 

Code is built upon previous work from Shady Rahayel (https://github.com/srahayel/SIR_mouse ), but modified to include the non-uniform clearance of aSyn from each brain region (defined by gene expression), to optimize parameters that best fit our output pathology maps, and to distinguish between aSyn accumulation and downstream atrophy (based on the methods from Zheng et al., 2019). 

## Main script: abm_optuna_general.py

## Full Workflow 

### Filter and package clearance genes into individual model inputs

1. **package_data.R**: provided by Yohan Yee, PhD, to parcellate the coronal and sagittal gene expression patterns for all genes from the Allen Mouse Brain Atlas' in-situ hybridization (ISH) data into the regions in the Allen Brain Ontology. Of particular interest are columns denoting the proportion of "valid" voxels containing an ISH measurement within each brain region, as well as the gene expression within each brain region normalized against whole-brain expression.
2. **filter_parcellated_GE.R**: find which regions have <20% valid voxels for each gene, and get rid of regions that are invalid across all genes (-2 regions). Additionally, discard any genes with >10/209 brain regions containing <20% valid voxels, indicating insufficient spatial coverage for a whole-brain analysis. For the remaining genes, fill in the "missing" regions with <20% valid voxels with the average gene expression across the whole brain. 
3. **pickle_abm_inputs.py**

### Run SIR agent-based model to find optimal parameters and simulate pathology 

4. **abm_optuna_general.py**: see above
5. **abm_clearance_genes.py**: allows you to simulate the "I" fraction/atrophy given a specific set of parameters for v, spread_rate, k1, k2, and retro. Output is the .pkl file containing the simulated neurodegeneration data for every timestep of the SIR model, but with a single set of "trial" parameters (no optuna). This is useful for simulating neurodegeneration in each region for the "optimal" parameterization determined above, and then correlating with the empirical pathology.

### Visualize simulated vs. empirical atrophy/pSyn pathology
6. **sir_conf_int.R**: come up with confidence intervals across 40 repeats of 200 Optuna trials for each baseline (without clearance), allowing us to determine the baseline correlation between simulated vs. empirical pathology without any clearance. This also plots the SIR model correlation after simulations using -t 1000, -t 3000, and -t 5000, allowing us to determine the number of time steps needed for the model to meaningfully converge.
7. **save_top_genes_and_venn.R**: look at a venn diagram of all genes that outperform the baselines from (6) for each epicentre/type of pathology
8.  **plot_avg_optuna.R**: plots smoothed correlation between simulated vs. empirical pathology across full range of parameters tested using Optuna, averaging across the top genes/gene ontology processes from the genes in (7). Need an input folder of .csv files generating using abm_optuna_general.py with the -x argument
9. **plot_sim_vs_actual_with_clearance.R**: plot the simulated vs. actual map of pathology, coloring each point by the relative expression of the clearance gene

## Miscellaneous other scripts

**plot_shady_ihc.R**: given templates/masks for the Allen Mouse Brain Atlas label file, as well as the regional average IHC-derived pSyn pathology from Rahayel et al., 2022, *Brain* with the same set of regions, plot brain maps of pSyn pathology to facilitate visual comparison with MRI-derived atrophies. 

**convert_shady_ihc_to_csv.R**: used to convert supplementary tables of IHC pathology from Rahayel et al., 2022, *Brain*.

## Dependencies
python 3.11.5 with all prerequisite packages installed
R with all prerequisite packages installed
qbatch is not necessary but highly recommended to parallelize job submission (especially when repeating abm_optuna_general.py with thousands of clearance genes and hundreds of trials)


## Results

See link to on-demand oral presentation and poster from ADPD 2025.
