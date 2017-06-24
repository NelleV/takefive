import os
from glob import glob
import argparse

import numpy as np
from sklearn.metrics import euclidean_distances
from iced.io import load_lengths


parser = argparse.ArgumentParser()
parser.add_argument("directory")
args = parser.parse_args()

filenames = glob(os.path.join(args.directory, "*NB0cst*_structure.txt"))
filenames.sort()

lengths = load_lengths(os.path.join(
    args.directory.replace("results", "data"),
    "raw.bed"))

distances = [[] for _ in lengths]

for filename in filenames:
    begin, end = 0, 0
    for l, length in enumerate(lengths):
        end += length
        X = np.loadtxt(filename)
        X = X[begin:end]
        mask = np.isnan(X[:, 0])
        X[mask] = 0
        dis = euclidean_distances(X)
        dis[mask] = np.nan
        dis[:, mask] = np.nan
        distances[l].append(dis[np.triu_indices(dis.shape[0], 1)][np.newaxis])
        begin = end

for l, distances_chr in enumerate(distances):
    distances_chr = np.concatenate(distances_chr)
    np.save(os.path.join(args.directory, "distances_chr%02d.npy" % l),
            distances_chr)
