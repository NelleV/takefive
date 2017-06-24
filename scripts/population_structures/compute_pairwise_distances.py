from __future__ import print_function

import os
from glob import glob
import argparse

import numpy as np
from sklearn.metrics import euclidean_distances


parser = argparse.ArgumentParser()
parser.add_argument("directory")
args = parser.parse_args()

filenames = glob(args.directory + "*_structure.txt")
filenames.sort()
filenames = filenames[:100]

distances = []
for filename in filenames:
    X = np.loadtxt(filename)
    mask = np.isnan(X[:, 0])
    X = X[np.invert(mask)]
    dis = euclidean_distances(X)
    distances.append(dis[np.triu_indices(dis.shape[0], 1)][np.newaxis])

distances = np.concatenate(distances)
outfile = args.directory.replace("structures", "results") + "_distances.npy"
print("Saving in %s" % outfile)

try:
    os.makedirs(os.path.dirname(outfile))
except OSError:
    pass

np.save(outfile, distances)
