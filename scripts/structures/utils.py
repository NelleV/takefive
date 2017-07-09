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

import numpy as np
from scipy import interpolate
from scipy.misc import comb


def bernstein_poly(i, n, t):
    """
     The Bernstein polynomial of n, i as a function of t
    """

    return comb(n, i) * (t**(n-i)) * (1 - t)**i


def bezier_curve(points, nTimes=100):
    """
       Given a set of control points, return the
       bezier curve defined by the control points.

       points should be a list of lists, or list of tuples
       such as [ [1,1],
                 [2,3],
                 [4,5], ..[Xn, Yn] ]
        nTimes is the number of time steps, defaults to 1000

        See http://processingjs.nihongoresources.com/bezierinfo/
    """

    nPoints = len(points[0])
    xPoints = np.array(points[0])
    yPoints = np.array(points[1])
    zPoints = np.array(points[2])

    f_z = interpolate.interp1d(
        np.arange(len(zPoints)), points[2])

    t = np.linspace(0.0, 1.0, nTimes)

    polynomial_array = np.array(
        [bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)])

    xvals = np.dot(xPoints, polynomial_array)
    yvals = np.dot(yPoints, polynomial_array)
    zvals = f_z(np.arange(len(xvals)))

    return xvals, yvals, zvals


def interpolate_chromosomes(X, lengths, kind="cubic", eps=1e-1, step=0.1):
    """
    Return a smoothed interpolation of the chromosomes

    Parameters
    ----------
    X : ndarray (n, 3)
        The 3D structure

    lengths : ndarray (L, )
        The lengths of the chromosomes. Note that the sum of the lengths
        should correspond to the length of the ndarray of the 3D structure.

    Returns
    -------
    smoothed_X : the smoothed 3D structure

    """
    smoothed_X = []
    mask = np.invert(np.isnan(X[:, 0]))

    begin, end = 0, 0
    for j, length in enumerate(lengths):
        end += length
        x = X[begin:end, 0]
        y = X[begin:end, 1]
        z = X[begin:end, 2]
        indx = mask[begin:end]

        if not len(x):
            break

        if not indx[0]:
            x[0] = x[indx][0]
            x[indx][0] = np.nan
            y[0] = y[indx][0]
            y[indx][0] = np.nan

            z[0] = z[indx][0]
            z[indx][0] = np.nan

        if not indx[-1]:
            x[-1] = x[indx][-1]
            x[indx][-1] = np.nan
            y[-1] = y[indx][-1]
            z[indx][-1] = np.nan

            z[-1] = z[indx][-1]
            z[indx][-1] = np.nan

        if kind is None:
            smoothed_X.append([x, y, z])
            begin = end
            continue

        indx = np.invert(np.isnan(x))
        m = np.arange(len(x))[indx]

        if kind == "bezier":
            m = np.arange(m.min(), len(x), step)
            f_x, f_y, f_z = bezier_curve(
                [x[indx], y[indx], z[indx]])
            smoothed_X.append([f_x,
                               f_y,
                               f_z])
        elif kind == "rbf":
            f_x = interpolate.Rbf(m, x[indx])
            f_y = interpolate.Rbf(m, y[indx])
            f_z = interpolate.Rbf(m, z[indx])

            m = np.arange(m.min(), len(x), 0.1)

            smoothed_X.append([f_x(np.arange(m.min(), m.max(), step)),
                               f_y(np.arange(m.min(), m.max(), step)),
                               f_z(np.arange(m.min(), m.max(), step))])

        else:
            f_x = interpolate.interp1d(m, x[indx],
                                       kind=kind)
            f_y = interpolate.interp1d(m, y[indx],
                                       kind=kind)
            f_z = interpolate.interp1d(m, z[indx],
                                       kind=kind)

            m = np.arange(m.min(), len(x), step)

            smoothed_X.append([f_x(np.arange(m.min(), m.max(), step)),
                               f_y(np.arange(m.min(), m.max(), step)),
                               f_z(np.arange(m.min(), m.max(), step))])

        begin = end
    return smoothed_X


def split_on_positions(chrm, chrm_num, gene_clusters, lengths,
                       resolution=10000):
    """
    """
    values = chrm_num * np.ones(len(chrm))
    for i, (chrom_cluster, begin, end) in enumerate(gene_clusters):
        if (chrom_cluster - 1) != chrm_num:
            continue
        length = lengths[chrm_num]
        to_length = float(len(chrm))
        res = float(length) / to_length
        end /= res * resolution
        begin /= res * resolution
        values[int(begin):int(end)] = len(lengths)
    return values
