import numpy as np
from scipy import sparse
import os
import argparse
from pastis import fastio
from pastis.optimization import poisson_structure
from pastis.optimization import mds
from utils import compute_wish_distances
import iced
from sklearn.metrics import euclidean_distances

"""
Launches the inference of the 3D model on .matrix/.bed files"""


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--seed", default=1, type=int)
parser.add_argument("--lengths", default=None, type=str)
args = parser.parse_args()

if args.filename.startswith("data"):
    filename = args.filename.replace("data", "results")
else:
    filename = args.filename

outname = filename.replace(
    ".matrix", "_PO_%02d_structure.txt" % (args.seed))
if os.path.exists(outname):
    import sys
    #sys.exit(0)

try:
    os.makedirs(os.path.dirname(filename))
except OSError:
    pass

try:
    lengths = fastio.load_lengths(args.filename.replace(".matrix", ".bed"))
except IOError:
    lengths = fastio.load_lengths(args.filename.replace(".matrix", "_abs.bed"))

counts = fastio.load_counts(args.filename, lengths=lengths)

# Remove the diagonal and remove 0 from matrix
counts.setdiag(0)
counts = counts.tocsr()
counts.eliminate_zeros()
counts = counts.tocoo()

print("Normalizing")

counts = counts.toarray()
counts = counts + counts.T
counts[np.arange(len(counts)), np.arange(len(counts))] = 0

counts = iced.filter.filter_low_counts(
    counts, lengths=lengths,
    sparsity=False, percentage=0.03)
normed, bias = iced.normalization.ICE_normalization(counts,
                                                    output_bias=True)

counts[np.isnan(counts)] = 0
normed[np.isnan(normed)] = 0
counts = sparse.coo_matrix(np.triu(counts))
normed = sparse.coo_matrix(np.triu(normed))
bias = bias.flatten()
bias[np.isnan(bias)] = 1
counts.setdiag(0)
counts = counts.tocsr()
counts.eliminate_zeros()
counts = counts.tocoo()
normed.setdiag(0)
normed = normed.tocsr()
normed.eliminate_zeros()
normed = normed.tocoo()

# Compute starting point
print("Initializing")
random_state = np.random.RandomState(args.seed)
wd = compute_wish_distances(normed)
X = mds.estimate_X(wd, random_state=random_state,
                   precompute_distances="precomputed")
X[np.isnan(X)] = 0

# PM2
outname = filename.replace(
    ".matrix", "_PO_%02d_structure.txt" % (args.seed, ))

if os.path.exists(outname):
    print("Already computed")
    import sys
    #sys.exit()

alpha = -3.
max_iter = 1000000

bias = bias.reshape(-1, 1)
counts = counts.toarray()
m  = (counts.sum(axis=0) + counts.sum(axis=1)) == 0
counts[m, :] = np.nan
counts[:, m] = np.nan

m, n = X.shape
d = euclidean_distances(X)
mask = (np.invert(np.tri(m, dtype=np.bool)) & np.invert(np.isnan(counts)))
beta = counts[mask].sum() / ((d[mask] ** alpha) * (bias * bias.T)[mask]).sum()
counts[np.isnan(counts)] = 0

counts = sparse.coo_matrix(counts)
X = poisson_structure.estimate_X(
    counts, alpha=-3., beta=beta, bias=bias,
    verbose=1,
    maxiter=max_iter,
    use_zero_entries=True,
    ini=X.flatten())

counts = np.array(counts.todense())
mask = (np.array(counts.sum(axis=0)) +
        np.array(counts.sum(axis=1)) == 0).flatten()
X[mask] = np.nan

np.savetxt(outname, X)
print("Finished", outname)
