#!/bin/bash

STATES=(rings trophozoites schizonts)
# FOLDERS=(HP1_KD HP1_tagged pv_sporozoite)
for STATE in ${STATES[@]}; do

  sbatch --export=STATE=${STATE} -v -p RM-shared --array 1-100 --mem 1500Mb \
    --job-name ${STATE}_VE \
    --time=4:00:00 cluster_scripts/slurm_ve.sh \
    -N 1 \
    --mail-type=ALL --mail-user=nelle@berkeley.edu;
done;


