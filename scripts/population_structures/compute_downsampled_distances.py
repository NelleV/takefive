import os
import sys
from glob import glob
import argparse

import numpy as np
from sklearn.metrics import euclidean_distances
from iced.utils import downsample_resolution
from iced.io import load_lengths


parser = argparse.ArgumentParser()
parser.add_argument("directory")
parser.add_argument("--lengths", "-l")
parser.add_argument("--factor", type=int, default=10)
args = parser.parse_args()

factor = args.factor
if args.lengths is not None:
    lengths = load_lengths(args.lengths)
else:
    lengths = load_lengths(os.path.join(
        args.directory.replace("results", "data"),
        "raw.bed"))

distances = []

filenames = glob(args.directory + "*_structure.txt")
filenames.sort()

distances = []

for i, filename in enumerate(filenames):
    sys.stdout.write(
        "\rAnalysing %0.2f %% files" % (100. * (i + 1) / len(filenames)))
    sys.stdout.flush()
    try:
        X = np.loadtxt(filename)
    except ValueError:
        print(i, filename)
        raise ValueError
    mask = np.isnan(X[:, 0])
    X[mask] = 0
    dis = euclidean_distances(X)
    dis[mask] = np.nan
    dis[:, mask] = np.nan
    dis, lengths_ = downsample_resolution(dis, lengths, factor=factor)
    distances.append(dis[np.triu_indices(dis.shape[0])][np.newaxis])

distances = np.concatenate(distances)
outfile = (
    args.directory.replace(
        "structures", "results") + "_distances_downsampled.npy")
print("Saving in %s" % outfile)

try:
    os.makedirs(os.path.dirname(outfile))
except OSError:
    pass

np.save(outfile, distances)
