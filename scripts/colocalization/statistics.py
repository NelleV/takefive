from scipy import spatial
from sklearn.metrics import euclidean_distances
import utils
import numpy as np


def compute_null_(chr_, X, lengths, nsamples=10000, shift=False,
                  method="mepd", thres=1):
    # Will be true except for same chromosomes
    m = spatial.distance.pdist(chr_[:, np.newaxis]).astype(bool)

    indices = np.zeros((len(chr_), nsamples))
    begin_chr = np.concatenate([[0], lengths.cumsum()])

    for i, l in enumerate(chr_):
        idx = np.random.choice(np.arange(lengths[l]), size=nsamples)
        indices[i] = idx + begin_chr[l]

    null = np.zeros(nsamples)
    indices = indices.astype(int).T
    for i, idx in enumerate(indices):
        x = X[idx]
        if method == "mepd":
            null[i] = np.median(spatial.distance.pdist(x)[m])
        else:
            null[i] = np.sum(spatial.distance.pdist(x)[m] < thres)
    return null


def compute_mepd(X, lengths, positions, chrom, return_null=False):
    assert len(positions) == len(chrom)
    positions = positions.copy()
    # Interpolate the missing values
    chromosomes = utils.interpolate_chromosomes(
        X, lengths, step=1,
        kind="rbf")

    chromosomes = [np.array(c) for c in chromosomes]
    chromosomes = np.concatenate(chromosomes, axis=1).T
    X = chromosomes

    # Create null distribution
    nsamples = 100000
    null = compute_null_(chrom, X, lengths, nsamples=nsamples, method="mepd")

    begin_chr = np.concatenate([[0], lengths.cumsum()])
    positions += begin_chr[chrom]

    m = spatial.distance.pdist(chrom[:, np.newaxis]).astype(bool)

    MEPD = np.median(spatial.distance.pdist(X[positions])[m])
    pval_l = 1. * (MEPD < null).sum() / nsamples
    pval_r = 1. * (MEPD > null).sum() / nsamples

    if return_null:
        return pval_l, pval_r, MEPD, null
    else:
        return pval_l, pval_r


def compute_witten_and_noble(X, lengths, positions, chrom, return_null=False,
                             threshold=0.10):
    assert len(positions) == len(chrom)
    positions = positions.copy()
    # Interpolate the missing values
    chromosomes = utils.interpolate_chromosomes(
        X, lengths, step=1,
        kind="rbf")

    chromosomes = [np.array(c) for c in chromosomes]
    chromosomes = np.concatenate(chromosomes, axis=1).T

    # There are more efficient ways to compute this.
    X[np.isnan(X)] = 0
    max_dis = euclidean_distances(X).max()
    X = chromosomes

    # Create null distribution
    nsamples = 100000
    null = compute_null_(
        chrom, X, lengths, nsamples=nsamples,
        method="witten", thres=(max_dis * threshold))

    begin_chr = np.concatenate([[0], lengths.cumsum()])
    positions += begin_chr[chrom]

    m = spatial.distance.pdist(chrom[:, np.newaxis]).astype(bool)

    MEPD = np.sum(
        spatial.distance.pdist(X[positions])[m] < max_dis * threshold)
    pval_l = 1. * (MEPD < null).sum() / nsamples
    pval_r = 1. * (MEPD > null).sum() / nsamples

    if return_null:
        return pval_l, pval_r, MEPD, null
    else:
        return pval_l, pval_r
