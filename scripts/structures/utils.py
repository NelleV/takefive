from scipy import sparse
import numpy as np


def compute_wish_distances(counts, alpha=-3.0, beta=1., use_zero_counts=True):
    counts = counts.toarray()
    m = (counts.sum(axis=0) + counts.sum(axis=1)) == 0
    wish_distances = counts.copy() / beta
    wish_distances[wish_distances != 0] **= 1. / alpha
    if use_zero_counts:
        wish_distances[counts == 0] = wish_distances[counts != 0].max()
    wish_distances[m] = 0
    wish_distances[:, m] = 0
    return sparse.coo_matrix(np.triu(wish_distances, 1))
