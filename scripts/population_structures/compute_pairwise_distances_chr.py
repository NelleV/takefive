import os
from glob import glob
import argparse

import numpy as np
from sklearn.metrics import euclidean_distances
from iced.io import load_lengths


parser = argparse.ArgumentParser()
parser.add_argument("directory")
parser.add_argument("--lengths")
args = parser.parse_args()

if args.lengths is not None:
    lengths = load_lengths(args.lengths)
else:
    lengths = load_lengths(os.path.join(
        args.directory.replace("results", "data"),
        "raw.bed"))

filenames = glob(args.directory + "*_structure.txt")
filenames.sort()
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
    outfile = (
        args.directory.replace(
            "structures", "results") + "_distances_chr%02d.npy" % (l+1))
    print("Saving in %s" % outfile)

    np.save(outfile,
            distances_chr)
