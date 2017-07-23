import argparse
import os
import numpy as np
from minorswing import fastio
import iced
import matplotlib.pyplot as plt
from matplotlib import colors


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--lengths")
parser.add_argument("--no-normalization", dest="no_normalization",
                    default=False, action="store_true")
args = parser.parse_args()
if args.lengths is not None:
    lengths = fastio.load_lengths(args.lengths)
else:
    try:
        lengths = fastio.load_lengths(args.filename.replace(".matrix", ".bed"))
    except IOError:
        try:
            lengths = fastio.load_lengths(
                args.filename.replace(".matrix", "_abs.bed"))
        except IOError:
            try:
                lengths = fastio.load_lengths(
                    args.filename.replace("_normalized.matrix", "_abs.bed"))
            except IOError:
                lengths = fastio.load_lengths(
                   args.filename.replace("_normalized.matrix", ".bed"))

counts = fastio.load_counts(args.filename, lengths=lengths)

centromeres = np.loadtxt("files/pf.cent")

AP2 = np.loadtxt(
    'files/AP2_all_positions.tab',
    usecols=(1, 2, 3), dtype=str)
AP2[:, 0] = [int(i[3:]) for i in AP2[:, 0]]
AP2 = AP2.astype(float)

if "25kb" in args.filename:
    resolution = 25000
elif "20000" in args.filename:
    resolution = 20000
else:
    resolution = 10000

counts = np.array(counts.todense())
counts = counts + counts.T
counts[np.arange(len(counts)), np.arange(len(counts))] = 0

if not args.no_normalization:
    counts = iced.filter.filter_low_counts(counts, lengths=lengths,
                                           sparsity=False, percentage=0.03)
    # counts = filter_low_mappability(counts, mappabilities=mappability,
    #                                 cutoff=0.25)
    counts, bias = iced.normalization.ICE_normalization(counts,
                                                        output_bias=True)

    counts[np.isnan(counts)] = 0

c = counts.copy().flatten()
c.sort()
vmax = c[int(0.99*len(c))]
c = c[c != 0]
vmin = c[int(0.05*len(c))]
fig, ax = plt.subplots(figsize=(12, 12))
m = ax.matshow(counts, origin="bottom", norm=colors.SymLogNorm(5),
               extent=(0, len(counts), 0, len(counts)),
               cmap="Blues",
               vmin=vmin,
               vmax=vmax)


if "pv" not in args.filename.lower():
    centromeres /= resolution
    centromeres = centromeres.mean(axis=1)
    centromeres[1:] += lengths[:-1].cumsum()
    [ax.axhline(i, linestyle="--", linewidth=1, color="#000000")
     for i in centromeres]
    [ax.axvline(i, linestyle="--", linewidth=1, color="#000000")
     for i in centromeres]
    lencum = np.concatenate([[0], lengths.cumsum()])
    ap2 = (AP2[:, 1:].mean(axis=1) / resolution +
           lencum[AP2[:, 0].astype(int) - 1])
    [ax.axhline(i, linestyle=":", linewidth=1, color="g") for i in ap2]
    [ax.axvline(i, linestyle=":", linewidth=1, color="g") for i in ap2]


ax.axhline(-0, color="#000000", linewidth=1, linestyle="-")
ax.axvline(-0, color="#000000", linewidth=1, linestyle="-")

cb = fig.colorbar(m)
cb.set_label("Contact counts")

ax.set_xlim((-0, len(counts) - 0))
ax.set_ylim((-0, len(counts) - 0))

if lengths is not None:
    [ax.axhline(i, linestyle="-", linewidth=1, color="#000000")
     for i in lengths.cumsum()]
    [ax.axvline(i, linestyle="-", linewidth=1, color="#000000")
     for i in lengths.cumsum()]

outname = args.filename.replace(
    "data", "images", 1).replace(".matrix", ".png")

try:
    os.makedirs(os.path.dirname(outname))
except OSError:
    pass

try:
    fig.savefig(outname, dpi=400)
except IOError:
    pass

centromeres = np.loadtxt("files/pf.cent")
centromeres = centromeres.mean(axis=1) / resolution
b, e = 0, 0
for i, l in enumerate(lengths):
    e += l
    scounts = counts[b:e, b:e]
    fig, ax = plt.subplots(figsize=(5, 5))
    m = ax.matshow(scounts, origin="bottom", norm=colors.SymLogNorm(5),
               extent=(0, l, 0, l),
               cmap="Blues")
    ax.axhline(centromeres[i], linestyle="--", linewidth=1, color="#000000")
    ax.axvline(centromeres[i], linestyle="--", linewidth=1, color="#000000")
    ax.set_xticks([])
    ax.set_yticks([])

    ap2 = AP2[:, 1:].mean(axis=1)[AP2[:, 0] - 1 == i] / resolution
    [ax.axhline(a, linestyle=":", linewidth=1, color="g") for a in ap2]
    [ax.axvline(a, linestyle=":", linewidth=1, color="g") for a in ap2]

    try:

        fig.savefig(
            outname.replace(
                ".pdf", "_%02d.pdf" % (i+1)).replace(
                    ".png", "_%02d.png" % (i+1)),
            dpi=400)
    except IOError:
        pass

    b = e
