import sys
import re
import numpy as np
import pandas as pd
import pickle
sys.path.insert(1, './SIR_mouse/model/')
from AgentBasedModel import AgentBasedModel
from scipy.stats import zscore, norm
from tqdm import tqdm
import argparse
import optuna
import scipy.stats

#Hipp injection (see find_peak_fit())
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--atrophy", type=str2bool, default=False, dest="atrophy",
                        nargs='?', help="Boolean for whether we want to simulate atrophy or stop at I fraction")
    parser.add_argument("-m", "--map_to_predict", dest="map",
                        nargs='?', help="Pathological map to compare to simulated atrophy")
    parser.add_argument("-x", "--suffix", dest="suffix", default=None,  nargs='?',
                        help="Suffix (typically used for multiple random repeats of same configuration)")
    parser.add_argument('-k', "--connectome", dest="connectome", nargs='?',
                        help="Connectome (pkl file w/ weights, distance, region size, sources, targets )")
    parser.add_argument(
        "-g", "--GEdir", default=None, dest="clearance_gene_dir",
         nargs='?', help="Clearance Gene Directory"
    )
    parser.add_argument(
        "-c", "--clearance", default=None, dest="clearance_gene",
        help="Specify the gene modulating clearance"
    )
    
    parser.add_argument(
        "-r", "--retro", default=False, dest="retro",
        type=str2bool, nargs='?', help="Retrograde spreading (default False) "
    )

    parser.add_argument(
        "-t", "--time", default=30000, type=int, nargs='?',
        dest="total_time", help="Total spreading time"
    )
    parser.add_argument(
        "-d", "--delta-t", default=0.1, type=float, nargs='?',
        dest="dt", help="Size of time increment"
    )
    parser.add_argument(
        "-e", "--epsilon", default=1e-5, type=float, nargs='?', #Shady's default was 1e-7 but this was waaaay too slow
        dest="eps", help="Convergence threshold"
    )
    
    parser.add_argument(
        "-S", "--seed", type=int, nargs='?',
        dest="seed", help="Simulated seeding site of misfolded alpha-syn"
        # injecting into the CP; CP = 35; CA1 = 24
    )
    parser.add_argument(
        "-y", "--num_trials", type=int, nargs='?', default=200,
        dest="num_trials", help="Number of optuna trials"
        # injecting into the CP; CP = 35; CA1 = 24
    )
    parser.add_argument(
        "-n", dest="study_name",
         nargs='?', help="Study name")
    args = parser.parse_args()
    return args

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def load_params(connectome):
 
    # load homogeneous values
    # choose the rate from 0.1 to 0.9
    #homorate = 0.5
    #syngene = np.transpose(np.full((1,213),homorate))
    with open(connectome, "rb") as f:
        params = pickle.load(f)
        # if retro is False, load data from 'params_nature.pickle'; anterograde spreading
        weights = params['weights']
        distance = params['distance']
        region_size = params['region_size']
        sources = params['sources']
        targets = params['targets']

    #Load Yohan's SNCA
    ge = pd.read_pickle('/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/Snca.pkl')

    syngene = ge.loc[sources, "Snca"] #can modify to be another gene if needed... 
    syngene = norm.cdf(zscore(syngene))

    return (
            weights, distance, np.append(region_size, region_size),
            sources, targets, np.append(syngene[:len(sources)], syngene[:len(sources)])
        )
         
def load_clearance(clearance_gene=None): 
    if clearance_gene is None:
        return 0.5 # default clearance rate
        # if clearance_gene is not specified (i.e., it is None), the function returns the default clearance rate of 0.5
    else:
        #CHANGE THIS: preprocessed .pkl file, need to COPY to both hemispheres and get rid of regions
        ge = pd.read_pickle(clearance_gene_dir+clearance_gene+'.pkl')
        epr=ge.loc[:,clearance_gene]
        #epr = ge.loc[sources, clearance_gene] 
        # #everything is in the right order - see SIR_mouse_exploration.ipynb on cicws        
        epr = np.append(epr, epr)
        # appends the data to itself

        return norm.cdf(zscore(epr.flatten()))
        # normalizes and transforms the data, and returns the resulting array as a cumulative distribution function

#copy in Steph's code to predict atrophy from I fraction using local degeneration + deafferentation
def simulate_pathology(atrophy, connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps, k1_atrophy=None, k2_atrophy=None):
    #taken from main method of abm_clearance_genes.py
    #inefficient to parse twice but meh we move

    weights, distance, region_size, sources, targets, syngene = load_params(connectome)
    clearance_rate = load_clearance(clearance_gene)

    abm = AgentBasedModel(
        weights=weights, distance=distance, region_size=region_size,
        sources=sources, targets=targets, dt=1
    )
    # reads input arguments passed to the script through the command line
        # such as the spread rate, the injection amount, and the total time. 

    abm.set_growth_process(growth_rate=syngene)
    abm.set_clearance_process(clearance_rate=clearance_rate)
    abm.set_spread_process(v=v)
    abm.update_spread_process(spread_scale=spread_rate) #CHANGE THIS: look into v scale
    #CHANGE THIS: look into (lack??) of update_growth_process, update_clearance_process, update_trans_process

    # calls several functions to load the model parameters, 
        # the clearance rate of proteins, and to set up the growth, clearance, and spreading processes of the ABM

    # growth process
    #print("Begin protein growth process....")
    for t in range(30000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        #abm.clearance_step_thresholded()
        abm.s_spread_step()
        #print("Number of regions where S hasn't converged:", np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum())
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    abm.dt = 0.1
    for t in range(500000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        #abm.clearance_step_thresholded()
        abm.s_spread_step()
        #print("Number of regions where S hasn't converged:", np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum())
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    abm.dt = dt
    #print("Third loop")
    for t in range(1000000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        #abm.clearance_step_thresholded()
        abm.s_spread_step()
        #print("Number of regions where S hasn't converged:", np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum())
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > eps, 1, 0).sum() == 0:
            break
    # runs the protein growth process in three stages
    # In each stage, the ABM model computes the protein growth, clearance, and spreading steps until the growth process stops, 
    # which is detected by comparing the changes in protein concentrations between time steps

    # spread process
    #print("Begin protein spreading process...")
    #print("Inject infectious proteins into {}...".format(abm.targets[seed]))
    abm.injection(seed=seed, amount=injection_amount)

    for t in tqdm(range(total_time)):
        try:
            abm.growth_step()
            abm.clearance_step()
            #abm.clearance_step_thresholded()
            abm.trans_step()
            #abm.trans_step_thresholded()
            abm.s_spread_step()
            abm.i_spread_step()
            abm.record_history_region()
            abm.record_history_edge()
        except:
            print("Numeric exception (see above); moving onto next Optuna trial")
            break
        
    infected_region_time=abm.i_region_history
    # Total number of proteins per regions per timestep
 
    if atrophy:
        normal_region_time=abm.s_region_history
        dt=abm.dt

        total_proteins_region_time = normal_region_time + infected_region_time

        # Ratio of infected over total proteins per region at every timestep
        infected_proteins_ratio_region_time = infected_region_time / total_proteins_region_time

        ## Normalized measure of connectivity strength for every target regions
        # step 1: Calculate the sum of connectivity strength for each source region
        sum_strength = np.sum(weights, axis=1)
        ratio_weights = weights / np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)
    
        k1 = k1_atrophy # 0.5
        k2 = k2_atrophy # 0.5

        # Create the new matrix with dimensions (426, 426)
        new_ratio_weights = np.zeros((len(targets), len(targets)))
        # Copy the last 213 rows to the first 213 rows of the new matrix
        new_ratio_weights[:len(sources), :] = ratio_weights[-len(sources):, :]
        # Copy the first 213 rows to the last 213 rows of the new matrix
        new_ratio_weights[len(sources):, :] = ratio_weights[:len(sources), :]
        # neuronal loss caused by lack of input from neighboring regions
        ratio_cum = np.dot(new_ratio_weights, (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
        ratio_cum = np.hstack([np.zeros((len(targets), 1)), ratio_cum[:, :-1]])
        ratio_cum = (k2 * ratio_cum) + (k1 * (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
        simulated = np.cumsum(ratio_cum, axis=1)
    else:
        simulated=infected_region_time.T #do not normalize by total number of proteins
    return simulated
    
def find_peak_fit(simulated, map):
    #modify from find_opt_params_clearance.r
    input_dir="./"
    #copied from above

    #Hipp injection
    pathology=pd.read_pickle(map)

    sim_pathology_data = pd.DataFrame(simulated)
    yohan_full = pd.read_pickle(input_dir+'yohan_source_full.pkl')
    ABA_region_names = yohan_full
    nregions = ABA_region_names.shape[0]
    ABA_region_names = np.vstack((ABA_region_names,ABA_region_names)).flatten()

    #process simulated atrophy data to get it in same form as Steph's DBM atrophies
    sim_pathology_data.index=ABA_region_names

    alternating_values = np.repeat(["right","left"], nregions)
    
    # Create unique row names based on first and last columns 
    sim_pathology_data.index = [alternating_values[i]+" "+list(sim_pathology_data.index)[i] for i in range(len(sim_pathology_data.index))]
    results=pathology
    common_row_names = np.intersect1d(results.index, sim_pathology_data.index)


    # Sort matrices based on common row names
    empirical_common = results.loc[common_row_names]
    simulated_common = sim_pathology_data.loc[common_row_names, :]
    nsteps=simulated_common.shape[1]
    correlations = np.zeros(nsteps)
    # Calculate Spearman correlation with each column in sim_pathology_data
    for i in range(nsteps):
        correlations[i] = scipy.stats.spearmanr(empirical_common.loc[:], simulated_common.loc[:, i]).correlation
    ##START HERE: MAKE SURE THIS WORKS WITH CHOPPED UP ATROPHY MAPS each has one column)
    peak_corr = np.max(correlations)
    peak_timestep = np.argmax(correlations)
    return (peak_corr, peak_timestep)

class Objective: ###EDIT HERE
    def __init__(self, atrophy, dt, eps, seed, total_time, clearance_gene, map, connectome, clearance_gene_dir=None): #connectome can either be retro/antero
        # Hold this implementation specific arguments as the fields of the class.
        self.atrophy = atrophy
        self.dt=dt
        self.eps = eps
        self.dt = dt
        self.seed = seed
        self.total_time = total_time
        self.clearance_gene = clearance_gene
        self.clearance_gene_dir = clearance_gene_dir
        self.map = map
        self.connectome = connectome

    def __call__(self, trial):
        # Calculate an objective value by using the extra arguments.
            v = trial.suggest_float("v",1e-3, 1000, log=True)
            spread_rate = trial.suggest_float("spread_rate",0.0001, 1)
            injection_amount = trial.suggest_float("injection_amount",1,100) #hypothesis: will need more aSyn injected to produce same result w/ increased clearance

            if self.atrophy: 
                k1_atrophy=trial.suggest_float("k1_atrophy",0.001,1)
                k2_atrophy=trial.suggest_float("k2_atrophy",0.001,1)
                simulated=simulate_pathology(self.atrophy, self.connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps, k1_atrophy, k2_atrophy)
            else:
                simulated=simulate_pathology(self.atrophy, connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps)
                
            peak_corr,peak_timestep=find_peak_fit(simulated, map)
            self.peak_timestep = peak_timestep
            trial.set_user_attr("peak_timestep", int(self.peak_timestep))
            return peak_corr
# Callback to store the best Objective instance
best_objective = {"instance": None}

def store_best_objective(study, trial):
    if study.best_trial == trial:
        best_objective["instance"] = trial.user_attrs["peak_timestep"]

###RUN OPTUNA###
if __name__ == "__main__":

    ###load all arguments before calling SIR functions
    args = parse_arguments()
    atrophy = args.atrophy #boolean indicating whether or not we want atrophy
    retro = args.retro
    dt = args.dt
    eps=args.eps
    seed = args.seed
    total_time = args.total_time 
    clearance_gene = args.clearance_gene
    clearance_gene_dir = args.clearance_gene_dir
    map = args.map
    connectome = args.connectome

    study_name=args.study_name
    suffix = args.suffix #if you want a .csv file output but no dashboard
    num_trials = args.num_trials
    
    if study_name is None:
        study = optuna.create_study(direction="maximize")
    else:
        study = optuna.create_study(study_name=study_name, storage="sqlite:///optuna_dbs/"+re.sub(" ", "_", study_name)+".sqlite3",direction="maximize")  
                                    # Create a new study with database --> can change when you aren't looking at every single gene lol
    #study=optuna.create_study(direction="maximize")
    study.optimize(Objective(atrophy, dt, eps, seed, total_time, clearance_gene, map, connectome, clearance_gene_dir), n_trials=num_trials, callbacks=[store_best_objective]) #CHANGED FOR HIPP RETRO (rly slow)

    # Retrieve the best Objective instance
    if best_objective["instance"] is not None:
        peak_timestep = best_objective["instance"]
    if suffix is None and study_name is None:
        pass
    else:
        df = study.trials_dataframe()
        if study_name is None:
            df.to_csv('./optuna_csvs/'+str(seed)+'/ret'+str(retro)+'_'+suffix+'.csv')
        else:
            df.to_csv('./optuna_csvs/'+str(seed)+'/'+study_name+'_ret'+str(retro)+'.csv')

        df.to_csv('./optuna_csvs/'+str(seed)+'/ret'+str(retro)+"_"+str(suffix)+'.csv')
    if atrophy:
        print(clearance_gene, study.best_trial.value, study.best_trial.params['v'], study.best_trial.params['spread_rate'],
                        study.best_trial.params['injection_amount'], study.best_trial.params['k1_atrophy'],
                        study.best_trial.params['k2_atrophy'], peak_timestep, sep=",") # Show the best value.
    else:
        print(clearance_gene, study.best_trial.value, study.best_trial.params['v'], study.best_trial.params['spread_rate'],
                        study.best_trial.params['injection_amount'], peak_timestep, sep=",") # Show the best value.

