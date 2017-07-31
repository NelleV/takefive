import numpy as np

from glob import glob
import argparse
from sklearn.metrics import euclidean_distances

import matplotlib.pyplot as plt



############################
# Without VRSM

filenames = glob("results/without_VRSM/trophozoites_plasmodium_*.txt.score")
filenames.sort()

without_vrsm_scores = []
for filename in filenames:
    without_vrsm_scores.append(np.loadtxt(filename))
without_vrsm_scores = np.array(without_vrsm_scores)
print("Mean scores without VRSM genes", without_vrsm_scores.mean())

filenames = glob("results/with_VRMS/trophozoites_plasmodium_*.txt.score")
filenames.sort()

with_vrsm_scores = []
for filename in filenames:
    with_vrsm_scores.append(np.loadtxt(filename))
with_vrsm_scores = np.array(with_vrsm_scores)
print("Mean scores with VRSM genes", with_vrsm_scores.mean())

