#!/usr/bin/env/python
import numpy as np
import pandas as pd
import pickle
import sys
import os

if __name__ == "__main__":

    mr = sys.argv[1]
    vv = sys.argv[2]
    input_dir = "/scratch/m/mchakrav/vnathan/sir_extended/analysis/"
    with open(input_dir+'params_nature_yohan_ccfv3.pkl', "rb") as f:
        params = pickle.load(f)
    sources = params['sources']

    mr_vv_string='mr'+str(int(mr))+'vv'+str(vv)
    full_ge = pd.read_csv(input_dir+'combined_rgn_ge_df_filt_'+mr_vv_string+'_colmean.csv', index_col=0)
    #add abbreviations from Steph's dataset
    steph_allen_rgn_names = pd.read_csv(input_dir+'ABAnewatlas_final_withacros.csv')
    rgn_name_to_acronym_dict = dict(zip(steph_allen_rgn_names['Structure'], steph_allen_rgn_names['acronym']))
    full_ge.index = [rgn_name_to_acronym_dict[i] for i in full_ge.index]

    df=full_ge
    output_dir = '/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/'+mr_vv_string+'/'
    os.mkdir(output_dir)

    # Iterate through each column in the DataFrame
    for column in df.columns:
        # Create a new DataFrame with only the current column
        new_df = df[[column]]
        new_df = new_df.loc[sources,:] #ensure regions are in the same order as params_nature_yohan_ccfv3.pkl before feeding into abm_clearance_genes.py - same as retro
        # Define the output CSV file name
        output_csv = f'{column}.csv'
        
        # Save the new DataFrame to the CSV file
        new_df.to_csv(output_dir+output_csv)
        
        #print(f'Saved {output_csv}')
ge_names_df = pd.DataFrame({"gene":df.columns})
ge_names_df.to_csv(output_dir+'gene_names_'+mr_vv_string+'.csv',index=False)