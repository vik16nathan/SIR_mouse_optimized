#!/bin/bash

clearance_gene_dir="/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/"
clearance_list="./ca1_ret_ihc_top80.csv"
mapfile -t list < <(tail -n +1 "${clearance_list}")

#clearance_list="top50_cl_cp_ant.csv"

#completed=($(awk -F"," '{print $1}' clearance_cp_ant_results_eps1e-5.txt))
mpi=24
t=1000

gene="Pfkfb2" 
for gene in "${list[@]}"; do
#for i in {1..40}; do
    #clearance_gene=$(sed -n "${i}p" "${clearance_gene_dir}${clearance_list}")
    #clearance_gene=$(sed -n "${i}p" "${clearance_list}")

    ##repeat results, no clearance
    #echo "python3 abm_optuna_general.py -a False -g "/project/m/mchakrav/vnathan/sir_extended/derivatives/yohan_ge_filt/mr10vv0.2/" -m "/project/m/mchakrav/vnathan/sir_extended/preprocessed/shady_ihc/Total Pathology_24 MPIHIP injection.pkl" -k params_nature_retro_yohan_ccfv3.pkl  -r True -t 1000 -d 0.1 -e 1e-5 -S 24" >> ca1_ret_ihc_24mpi_t1000.csv
    #clearance gene
    echo "python3 abm_optuna_general.py -a False -c "${gene}" -g "${clearance_gene_dir}"  -m \"/project/m/mchakrav/vnathan/sir_extended/preprocessed/shady_ihc/Total Pathology_24 MPIHIP injection.pkl\" -k "params_nature_retro_yohan_ccfv3.pkl"  -r True -t $t -d 0.1 -e 1e-5 -S 24 -x "${gene}_ihc_24mpi" >> ca1_ret_ihc_top80_results.txt" >> ca1_ret_ihc_top80_commands.txt
done
   
