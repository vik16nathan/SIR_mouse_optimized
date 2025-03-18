# -*- coding: utf-8 -*-
"""This script runs the agent-based Susceptible-Infected-Removed (SIR) Model.
Authors:
    Ying-Qiu Zheng, Shady Rahayel
    
For running the model, run:
    
    python abm.py --retro True --speed 10 --spreading-rate 0.01 --time 30000
        --delta-t 0.1 --seed -1 --seed-amount 1
        
--retro True specifies a retrograde spreading
--speed is the spreading speed of agents in edges
--spreading-rate is the probability of staying inside a region
--time is the spreading time of agents
--delta-t is the size of timesteps
--seed is an integer that refers to the list of regions listed
alphabetically from the Allen Mouse Brain Atlas (see params_nature_retro.pickle)
CP = 35, ACB = 3, and CA1 = 24
--seed-amount is the initial injected amount of infected agents
        
This generates arrays containing the number of normal and infected agents
at each iteration for every region of the Allen Mouse Brain Atlas.
The distribution of normal agents can be found in .s_region_history
The distribution of infected agents can be found in .i_region_history

"""

import sys
import numpy as np
import pandas as pd
import pickle
sys.path.insert(1, './SIR_mouse/model/')
from AgentBasedModel import AgentBasedModel
from scipy.stats import zscore, norm
from tqdm import tqdm
import argparse
import time


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-g", "--GEdir", dest="clearance_gene_dir",
         nargs='?', help="Clearance Gene Directory"
    )

    parser.add_argument(
        "-o", "--output-dir", dest="output_dir",
        nargs='?', help="Result Output Directory"
    )
    
    parser.add_argument(
        "-r", "--retro", default=True, dest="retro",
        type=str2bool, nargs='?', help="Retrograde spreading (True) "
    )
    parser.add_argument(
        "-v", "--speed", default=10, dest="v", nargs='?',
        help="Spreading speed", type=float
    )
    parser.add_argument(
        "-s", "--spreading-rate", default=0.01, type=float,
        nargs='?', help="Spreading rate", dest="spread_rate"
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
        "-S", "--seed", default=35, type=int, nargs='?',
        dest="seed", help="Simulated seeding site of misfolded alpha-syn"
        # injecting into the CP; CP = 35; CA1 = 24
    )
    parser.add_argument(
        "-a", "--seed-amount", default="1", type=float,
        dest="injection_amount", help="Seeding amount of misfolded alpha-syn"
    )
    parser.add_argument(
        "-c", "--clearance", default=None, dest="clearance_gene",
        help="Specify the gene modulating clearance"
    )
    parser.add_argument(
        "-k1", "--k1", default=None, type=float, dest="k1_atrophy", nargs='?',
        help="k1_atrophy"
    )
    parser.add_argument(
        "-k2", "--k2", default=None, type=float, dest="k2_atrophy", nargs='?',
        help="k2_atrophy"
    )
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


def load_params(retro=False):
 
    # load homogeneous values
    # choose the rate from 0.1 to 0.9
    #homorate = 0.5
    #syngene = np.transpose(np.full((1,213),homorate))
    if retro is False:
        with open('./params_nature_yohan_ccfv3.pkl', "rb") as f:
            params = pickle.load(f)
            # if retro is False, load data from 'params_nature.pickle'; anterograde spreading
        #print(params)
        #FILTER TO GET RID OF INVALID REGIONS!
        weights = params['weights']
        distance = params['distance']
        region_size = params['region_size']
        sources = params['sources']
        targets = params['targets']

    elif retro is True: #IGNORE FOR NOW... 
        with open('./params_nature_retro_yohan_ccfv3.pkl', "rb") as f:
            params = pickle.load(f)
            # If retro is True, load data from 'params_nature_retro.pickle'; retrograde spreading
        
        weights = params['weights']
        distance = params['distance']
        region_size = params['region_size']
        sources = params['sources']
        targets = params['targets']

    #Load Yohan's SNCA
    ge = pd.read_csv('/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/Snca.csv',index_col=0)
    
    syngene = ge.loc[sources, "Snca"]
    syngene = norm.cdf(zscore(syngene))
    #Load Shady's SNCA
    #with open('/data/chamal/projects/natvik/sir_extended/analysis/shady_snca_rgn_filt.pkl', 'rb') as f:
    #    epr = pickle.load(f)
    #    syngene = epr['syngene']
    return (
            weights, distance, np.append(region_size, region_size),
            sources, targets, np.append(syngene[:len(sources)], syngene[:len(sources)])
        )
         
def load_clearance(clearance_gene=None): 
    if clearance_gene is None:
        return 0.5 # default clearance rate
        # if clearance_gene is not specified (i.e., it is None), the function returns the default clearance rate of 0.5
    else:
        #CHANGE THIS: preprocessed .csv file, need to COPY to both hemispheres and get rid of regions
        ge = pd.read_csv(clearance_gene_dir+clearance_gene+'.csv',index_col=0)
        epr=ge.loc[:,clearance_gene]
        #epr = ge.loc[sources, clearance_gene] 
        # #everything is in the right order - see SIR_mouse_exploration.ipynb on cicws        
        epr = np.append(epr, epr)
        # appends the data to itself

        return norm.cdf(zscore(epr.flatten()))
        # normalizes and transforms the data, and returns the resulting array as a cumulative distribution function


if __name__ == "__main__":
    # run ABM
    # read arguments
    args = parse_arguments()

    retro = args.retro
    #injection_site = args.injection_site
    v = args.v
    spread_rate = args.spread_rate
    dt = args.dt
    seed = args.seed
    injection_amount = args.injection_amount
    total_time = args.total_time
    clearance_gene = args.clearance_gene
    k1_atrophy = args.k1_atrophy
    k2_atrophy = args.k2_atrophy


    output_dir=args.output_dir
    clearance_gene_dir=args.clearance_gene_dir

    weights, distance, region_size, sources, targets, syngene = load_params(retro=retro)
    clearance_rate = load_clearance(clearance_gene)
    #with open('snca_norm.pickle', 'wb') as f:
    #    pickle.dump(syngene, f)
    #clearance_rate = load_clearance(clearance_gene) 
    #CHANGE THIS: add a filter step to get rid of regions that aren't represented in the clearance gene expression file
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
    print("Begin protein growth process....")
    start_time = time.time()
    for t in range(30000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > 1e-5, 1, 0).sum() == 0:
            break
    abm.dt = 0.1
    for t in range(500000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > 1e-5, 1, 0).sum() == 0:
            break
    abm.dt = dt
    for t in range(1000000000):
        prev = np.copy(abm.s_region)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.where(np.abs(prev - abm.s_region) / abm.s_region > 1e-5, 1, 0).sum() == 0:
            break
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Protein growth time: {elapsed_time:.4f} seconds")
    # runs the protein growth process in three stages
    # In each stage, the ABM model computes the protein growth, clearance, and spreading steps until the growth process stops, 
    # which is detected by comparing the changes in protein concentrations between time steps

    # spread process
    print("Begin protein spreading process...")
    start_time = time.time()
    print("Inject infectious proteins into {}...".format(abm.targets[seed]))
    abm.injection(seed=seed, amount=injection_amount)

    for t in tqdm(range(total_time)):
        abm.growth_step()
        abm.clearance_step()
        abm.trans_step()
        abm.s_spread_step()
        abm.i_spread_step()
        abm.record_history_region()
        abm.record_history_edge()
    
    # runs the protein spreading process for the specified total time. 
    # injects infectious proteins into the ABM's target region, 
    # simulates the protein growth, clearance, and spreading steps, 
    # and records the resulting protein concentrations at each time step. 

    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Protein spreading time: {elapsed_time:.4f} seconds")
    # saves the results of the ABM simulation in a pickle file

    # generates arrays containing the number of normal and infected agents at each iteration - LOOK AT ZHENG PAPER
    # for every region of the Allen Mouse Brain Atlas
        # The distribution of normal agents can be found in .s_region_history 
        # The distribution of infected agents can be found in .i_region_history

        # i_region_history = the number of infected agents in every region for every time step), 
        # s_region_history = the number of susceptible agents in every region for every time step)
        # i_edge_history = the number of infected agents in every edge for every time step), 
        # s_edge_history = the number of susceptibles agents in every edge for every time step)


    start_time = time.time()
    infected_region_time=abm.i_region_history
    #infected_edge=results.i_edge_history
    normal_region_time=abm.s_region_history
    #normal_edge=results.s_edge_history
    dt=abm.dt
    # size of timesteps

    # Total number of proteins per regions per timestep
    total_proteins_region_time = normal_region_time + infected_region_time

    # Ratio of infected over total proteins per region at every timestep
    infected_proteins_ratio_region_time = infected_region_time / total_proteins_region_time
    #ratio(Rmis_all<1) = 0; % remove possible NaNs... from MATLAB code

    if k1_atrophy is None and k2_atrophy is None:
        simulated=infected_region_time.T
        abm_filename = "abm_spread_v.{0}.spread_rate.{1}.dt.{2}.seed.{3}.injection_amount.{4}.clearance_gene.{5}." \
                .format(v, spread_rate, dt, seed, injection_amount, clearance_gene)
        np.savetxt(output_dir+abm_filename + '.csv', simulated, delimiter=',', fmt='%s')
        sys.exit(0)


    abm_filename = "abm_spread_v.{0}.spread_rate.{1}.dt.{2}.seed.{3}.injection_amount.{4}.clearance_gene.{5}.k1.{6}.k2.{7}" \
        .format(v, spread_rate, dt, seed, injection_amount, clearance_gene, k1_atrophy, k2_atrophy)

    if retro is True:
        abm_filename = "retro_" + abm_filename

    abm_filename_pickle = abm_filename + ".pickle"

    with open(output_dir+abm_filename_pickle, "wb") as f:
        pickle.dump(abm, f)

    ## Normalized measure of connectivity strength for every target regions
    # step 1: Calculate the sum of connectivity strength for each source region
    sum_strength = np.sum(weights, axis=1)
        # weights is defined as a matrix with dimensions N_regions by N_regions; 213 x 426 
        # np.sum(weights, axis=1); sum across every row to get total connectivity strength across all injection sites (target regions) for each region	

    # step 2: Calculate the ratio of weights to the sum of strength for each source region
    ratio_weights = weights / np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)
        # np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)
        # takes the sum of the connectivity strength for each source regions across every target region repeated N_target_regions times

    # # Print the resulting ratio weights
    # print(ratio_weights.shape) (213,426)
    #         # from this ratio, the resulting weights matrix provides a normalized measure of connectivity strength 
    #         # that takes into account the relative contribution of each source region to the overall connectivity with the target regions.

    k1 = k1_atrophy # 0.5
    #k2 = 1 - k1
    k2 = k2_atrophy # 0.5

    # input weights of deafferentation (scaled by structural connectivity)

    # Create the new matrix with dimensions (426, 426)
    new_ratio_weights = np.zeros((len(targets), len(targets)))
    # Copy the last 213 rows to the first 213 rows of the new matrix
    new_ratio_weights[:len(sources), :] = ratio_weights[-len(sources):, :]
    # Copy the first 213 rows to the last 213 rows of the new matrix
    new_ratio_weights[len(sources):, :] = ratio_weights[:len(sources), :]
    # print(new_ratio_weights.shape) #(426, 426)

    # neuronal loss caused by lack of input from neighboring regions
    ratio_cum = np.dot(new_ratio_weights, (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
            # apply a decay or scaling factor to the weights, modelling time-dependent changes or 
            # attenuation in connectivity strengths as a result of atrophy/deaff for example
    # print(ratio_cum.shape) #(426, 30000)

    # simulate a backward shift of one time step
    ratio_cum = np.hstack([np.zeros((len(targets), 1)), ratio_cum[:, :-1]])
            # shifts the ratio_cum matrix one time step back by 
                    # appending a column of zeros at the beginning and removing the last column of ratio_cum. 
            # The result is that each column in ratio_cum is shifted one position to the right, 
                    # with a column of zeros added at the leftmost position.
    # print(ratio_cum.shape) #(426, 30000)
    
    ratio_cum = (k2 * ratio_cum) + (k1 * (1 - np.exp(-infected_proteins_ratio_region_time.T * dt)))
    # update the ratio_cum matrix by applying a weighted combination of the shifted ratio_cum and a term involving k1, ratio, and dt. 
        # k2 * ratio_cum: This term scales the shifted ratio_cum matrix by a factor k2.
        # k1 * (1 - exp(-ratio * dt)): calculates the exponential decay factor 1 - exp(-ratio * dt) and scales it by k1.
        # The resulting ratio_cum matrix is obtained by adding these two terms together.
    # print(ratio_cum.shape) #(426, 30000)

    # add all the increments across each dt
    simulated_atrophy = np.cumsum(ratio_cum, axis=1)
        # np.cumsum(ratio_cum, axis=1) computes the cumulative sum along the second axis (axis=1) of the ratio_cum array. 
        # The resulting array simulated_atrophy will have the same shape as ratio_cum, and 
        # each element will contain the cumulative sum of the corresponding elements in ratio_cum.
    # print(simulated_atrophy) 
    # print(simulated_atrophy.shape) #(426, 30000)

    # Save variables to a file
    with open(output_dir+'saved_atrophy_' + abm_filename_pickle, 'wb') as f:
        pickle.dump((simulated_atrophy, infected_proteins_ratio_region_time, targets, dt), f)

    #save to a csv file

    #import rpy2.robjects as robjects
    np.savetxt(output_dir+abm_filename + '.csv', simulated_atrophy, delimiter=',', fmt='%s')
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Regional atrophy map post time: {elapsed_time:.4f} seconds")
