#! /bin/sh

source activate IMP
cd $SLURM_SUBMIT_DIR

python plasmo_landmark_without_VRSM.py ${SLURM_ARRAY_TASK_ID} --state ${STATE}
