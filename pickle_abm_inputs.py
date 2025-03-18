import numpy as np
import pandas as pd
import os
import glob
import multiprocessing

if __name__ == "__main__":
    #Pickle all output atrophy maps
    #csv_directory ='/project/m/mchakrav/vnathan/sir_extended/derivatives/'
    #pickle_directory ='/project/m/mchakrav/vnathan/sir_extended/derivatives/'
    #file_list_to_pickle = [
    #    "hipp_rgn_t_stats_hemiMsPff_hemiPBS.csv",
    #    "rgn_t_stats_full_hemiMsPff_hemiPBS.csv",
    #    "rgn_t_stats_full_hemiMsPff_wtPBS.csv"
    #]

    #pickle Shady's IHC data
    csv_directory = '/project/m/mchakrav/vnathan/sir_extended/derivatives/merfish/cell_class_counts/'
    pickle_directory = '/project/m/mchakrav/vnathan/sir_extended/derivatives/merfish/cell_class_counts/'
   # file_list_to_pickle = [
   #    'HIP injection.csv',
   #    'CP injection.csv',
   #    'ACB injection.csv'
   # ]

    #pickle gene expression
    #csv_directory = "/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/"
    #pickle_directory = "/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/"
    file_list_to_pickle = glob.glob(csv_directory + '*.csv')
    for file in file_list_to_pickle:
        #df = pd.read_csv(csv_directory + file, index_col=0)
        #for gene expression
        df = pd.read_csv(file, index_col=0)
        base = os.path.basename(file)
        pickle_filename = base.replace('.csv', '.pkl')
        pickle_path = os.path.join(pickle_directory, pickle_filename)
        df.to_pickle(pickle_path)
        # Create a pickle file name based on the CSV file name - UNCOMMENT FOR EVERYTHING THAT ISN'T GE
        #for column in df.columns:

            #base = column+os.path.basename(file)
            #pickle_filename = base.replace('.csv', '.pkl')
            #pickle_path = os.path.join(pickle_directory, pickle_filename)
            
            # Save the DataFrame as a pickle file
            #df[column].to_pickle(pickle_path)


    ##Pickle Yohan file names
    #df = pd.read_csv("yohan_source_full.csv", index_col=0)
    #df.to_pickle("yohan_source_full.pkl")