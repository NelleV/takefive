import os
import numpy as np
from sklearn.metrics import euclidean_distances
from mayavi import mlab
from minorswing import fastio
import argparse

from utils import interpolate_chromosomes
from minorswing.validation import realign_structures


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--lengths", default=None, type=str)
parser.add_argument("--counts", default=None, type=str)
parser.add_argument("--resolution", default=10000, type=int)
parser.add_argument("--outfile", "-o", default=None, type=str)
parser.add_argument("--rotate", "-r", default=False, action="store_true")
parser.add_argument("--nucleus-size", "-n", dest="nucleus_size", default=1000,
                    type=float)
parser.add_argument("--show", "-s", default=False, action="store_true")
parser.add_argument("--no-smooth", dest="smooth", default=True,
                    action="store_false")
parser.add_argument("--flipped", "-f", default=False, action="store_true")
args = parser.parse_args()

if "Pv" in args.filename:
    gene_clusters = np.loadtxt("files/pv_loci.txt", dtype=int)
else:
    gene_clusters = np.loadtxt(
        'files/VRSSMclusterIDs.txt',
        skiprows=1,
        usecols=(0, 2, 3), dtype=int)

if args.rotate and "Pv" not in args.filename:
    X_ = np.loadtxt(
        "results/gametocyte_structure.txt")

resolution = args.resolution

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

centromeres = np.loadtxt("files/pf.cent")

mlab.figure(bgcolor=(1, 1, 1))
mlab.clf()

f = mlab.gcf()
f.scene.light_manager.light_mode = 'vtk'

scale = 5
X = np.loadtxt(args.filename) * scale
if len(X.shape) != 2:
    X = X.reshape(-1, 3)

if "Pv" not in args.filename:
    # Realign:
    if args.rotate:
        mask = np.invert(np.isnan(X[:, 0]) | np.isnan(X_[:, 0]))
        X, _ = realign_structures(
            X_,
            X, rescale=True)
    else:
        mask = np.invert(np.isnan(X[:, 0]))
else:
    mask = np.invert(np.isnan(X[:, 0]))

# Rescale the structure to the proper size
# Divide nucleus size to look better
X *= args.nucleus_size / 100 / euclidean_distances(X[mask]).max()

if args.smooth:
    smoothed_X = interpolate_chromosomes(X, lengths, kind="rbf")
else:
    smoothed_X = interpolate_chromosomes(X, lengths, kind=None)

begin, end = 0, 0
for j, (x, y, z) in enumerate(smoothed_X):
    if args.flipped:
        z *= -1
    m = len(x)
    end += m
    length = lengths[j]
    cent = centromeres[j]

    mlab.plot3d(x, y, z,
                j * np.ones(x.shape), colormap="Spectral",
                vmin=0, vmax=14)
    shift = 0
    linewidth = 0.06
    phi, theta = np.mgrid[0:np.pi:101j, 0:2 * np.pi:101j]
    xp = 2 * linewidth * np.sin(phi) * np.cos(theta) + x[0]
    yp = 2 * linewidth * np.sin(phi) * np.sin(theta) + y[0]
    zp = 2 * linewidth * np.cos(phi) + z[0]
    mlab.mesh(xp + shift, yp, zp, opacity=1, color=(1, 1, 1))

    xp = 2 * linewidth * np.sin(phi) * np.cos(theta) + x[-1]
    yp = 2 * linewidth * np.sin(phi) * np.sin(theta) + y[-1]
    zp = 2 * linewidth * np.cos(phi) + z[-1]
    mlab.mesh(xp + shift, yp, zp, opacity=1, color=(1, 1, 1))

    for c, b, e in gene_clusters:
        if (c - 1) != j:
            continue
        else:
            pos = (b + e) / 2
            pos /= float(length) * resolution / m
            pos = pos.astype(int)
            xp = 2 * linewidth * np.sin(phi) * np.cos(theta) + x[pos]
            yp = 2 * linewidth * np.sin(phi) * np.sin(theta) + y[pos]
            zp = 2 * linewidth * np.cos(phi) + z[pos]
            mlab.mesh(xp + shift, yp, zp, opacity=1,
                      color=(0.07, 0.68, 0.37))

    if "Pv" in args.filename:
        begin = end
        continue

    ind = int(cent.mean() / resolution / length * m)
    xp = 2 * linewidth * np.sin(phi) * np.cos(theta) + x[ind]
    yp = 2 * linewidth * np.sin(phi) * np.sin(theta) + y[ind]
    zp = 2 * linewidth * np.cos(phi) + z[ind]
    mlab.mesh(xp + shift, yp, zp, opacity=1, color=(0.67, 0.77, 0.93))

    begin = end


if args.outfile is None:
    outfile = args.filename.replace(
        "results", "images").replace(".txt", ".png")
else:
    outfile = args.outfile

try:
    os.makedirs(os.path.dirname(outfile))
except OSError:
    pass
mlab.savefig(outfile, magnification=4)
if args.show:
    mlab.show()
