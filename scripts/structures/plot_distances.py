import os
import argparse

import numpy as np
from sklearn.metrics import euclidean_distances
import matplotlib.pyplot as plt

from minorswing import fastio


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--lengths", default=None, type=str)
parser.add_argument("--counts", default=None, type=str)

args = parser.parse_args()

if args.lengths is not None:
    if args.lengths.endswith(".npy"):
        lengths = np.load(args.lengths)
    else:
        lengths = fastio.load_lengths(args.lengths)
else:
    lengths = fastio.load_lengths(
        os.path.join(os.path.dirname(args.filename).replace(
            "results",
            "data"), os.path.basename(args.filename)[:13] + ".bed"))

if args.counts is not None:
    counts = fastio.load_counts(args.counts)
else:
    counts = None

centromeres = np.loadtxt("files/pf.cent")
resolution = 10000

X = np.loadtxt(args.filename)
if len(X.shape) == 1:
    X = X.reshape(-1, 3)

if counts is None:
    mask = np.isnan(X[:, 0])
else:
    mask = np.array((counts.sum(axis=0) == 0)).flatten()

X[np.isnan(X)] = 0
dis = euclidean_distances(X)
dis[mask] = np.nan
dis[:, mask] = np.nan

fig, ax = plt.subplots(figsize=(8, 8))
ax.matshow(dis, cmap="RdBu", shape=(0, len(dis), 0, len(dis)))
ax.set_xlim((-0, len(dis)))
ax.set_ylim((-0, len(dis)))
if lengths is not None:
    [ax.axhline(i, linestyle="--", linewidth=1, color="#000000")
     for i in lengths.cumsum()]
    [ax.axvline(i, linestyle="--", linewidth=1, color="#000000")
     for i in lengths.cumsum()]

if "pv" not in args.filename.lower():
    centromeres /= resolution
    centromeres = centromeres.mean(axis=1)
    centromeres[1:] += lengths[:-1].cumsum()
    [ax.axhline(i, linestyle=":", linewidth=1, color="#000000")
     for i in centromeres]
    [ax.axvline(i, linestyle=":", linewidth=1, color="#000000")
     for i in centromeres]

ax.axhline(-0.5, color="#000000", linewidth=1, linestyle="--")
ax.axvline(-0.5, color="#000000", linewidth=1, linestyle="--")

outfile = args.filename.replace("results", "images").replace(
    "structure.txt",
    "distances.png")
try:
    os.makedirs(os.path.dirname(outfile))
except OSError:
    pass

fig.savefig(outfile)
