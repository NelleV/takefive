#! /bin/sh

source activate minorswing
cd $SLURM_SUBMIT_DIR

filename=`data/ay2013/${DATASET}_10000_raw.matrix`
python infer_structures_nb.py  $filename \
  --seed ${SLURM_ARRAY_TASK_ID}
