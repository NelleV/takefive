#!/bin/bash

FOLDERS=(rings trophozoites schizonts)
# FOLDERS=(HP1_KD HP1_tagged pv_sporozoite)
for FOLDER in ${FOLDERS[@]}; do
#  sbatch --export=DATASET=${FOLDER} -v -p RM-shared --array 1-100 --mem 1500Mb \
#    --job-name ${FOLDER}_MDS \
#    --time=48:00:00 cluster_scripts/slurm_mds.sh \
#    -N 1 \
#    --mail-type=ALL --mail-user=nelle@berkeley.edu;
#
#  sbatch --export=DATASET=${FOLDER} -v -p RM-shared --array 1-100 --mem 1500Mb \
#    --job-name ${FOLDER}_NB0 \
#    --time=48:00:00 cluster_scripts/slurm_nb0.sh \
#    -N 1 \
#    --mail-type=ALL --mail-user=nelle@berkeley.edu;

  sbatch --export=DATASET=${FOLDER} -v -p RM-shared --array 101-1000 --mem 1500Mb \
    --job-name ${FOLDER}_PO \
    --time=48:00:00 cluster_scripts/slurm_po.sh \
    -N 1 \
    --mail-type=ALL --mail-user=nelle@berkeley.edu;


done;


