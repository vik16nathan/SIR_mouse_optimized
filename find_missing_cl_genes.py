import numpy as np
import pandas as pd

if __name__ == "__main__":
    
    # Load the clearance gene list
    clearance_genes = pd.read_csv('/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/gene_names_mr10vv0.2.csv')['gene']
    
    # Define file paths for command files and corresponding output CSVs
    command_filepaths = [
       'hipp_ant_ihc_clearance.txt', 
       'hipp_ret_ihc_clearance.txt'
    ]
    output_csvs = [
        'ca1_ant_ihc_clearance_24mpi_t1000.csv', 
        'ca1_ret_ihc_clearance_24mpi_t1000.csv'
    ]

    # Find which clearance genes in the overall list aren't in each output CSV
    for i in range(len(output_csvs)):
        out_csv = output_csvs[i]
        #print(out_csv)
        command_filepath = command_filepaths[i]
        
        # Read .txt file and store each line as an element in a list
        with open(command_filepath, 'r') as f:
            qbatch_lines = f.readlines()
        
        # Load genes present in the current output CSV
        complete_genes = pd.read_csv(out_csv, header=None)[0]
        # Iterate over each gene in the clearance gene list
        for gene in clearance_genes:
            if gene not in complete_genes.values:
                # Find and print the line from qbatch_lines containing the missing gene
                for line in qbatch_lines:
                    if gene in line:
                        print(line.strip())
                        break
