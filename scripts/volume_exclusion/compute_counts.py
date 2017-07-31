import numpy as np
from glob import glob
from sklearn.metrics import euclidean_distances
from joblib import Memory
from minorswing import fastio
from iced import utils
import iced
from utils import get_expected, get_mapping

from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt

mem = Memory(".joblib")

filenames = glob("results/without_VRSM/trophozoites_plasmodium_*.txt")
filenames.sort()

chr_seq = {"chr01": 640851, "chr02": 947102, "chr03": 1067971,
           "chr04": 1200490,
           "chr05": 1343557, "chr06": 1418242, "chr07": 1445207,
           "chr08": 1472805,
           "chr09": 1541735, "chr010": 1687656, "chr011": 2038340,
           "chr012": 2271494,
           "chr013": 2925236, "chr014": 3291936}

chr_id = [
    int(chr_[3:]) for chr_ in chr_seq.keys()
    for i in range(np.ceil(chr_seq[chr_] / 3200.).astype(int))]
lengths = [int(chr_seq["chr0%d" % int(d)]) for d in range(1, 10)] + \
          [int(chr_seq["chr0%d" % int(d)]) for d in range(10, 15)]
lengths = np.round(np.array(lengths) / 9600.).astype(int)


def build_counts(thres=45):
    # Only displaying chromosome 7
    counts = None

    for i, filename in enumerate(filenames):
        try:
            X = np.loadtxt(filename)
        except ValueError:
            print("Cannot load %s" % filename)
            continue
        X[:, 0] = chr_id
        X = X[X[:, 0] == 7, 1:]
        dis = euclidean_distances(X)
        if counts is not None:
            counts += dis <= thres
        else:
            counts = (dis <= thres).astype(int)
    return counts


counts = mem.cache(build_counts)(thres=75)
counts = counts.astype(float)
n = counts.shape[0]
idx = np.arange(n)
counts[idx[:-1], idx[1:]] = 0
counts[idx[1:], idx[:-1]] = 0

counts[idx[:-2], idx[2:]] = 0
counts[idx[2:], idx[:-2]] = 0

counts[idx[:-3], idx[3:]] = 0
counts[idx[3:], idx[:-3]] = 0
counts[idx[:-4], idx[4:]] = 0
counts[idx[4:], idx[:-4]] = 0
counts[idx[:-5], idx[5:]] = 0
counts[idx[5:], idx[:-5]] = 0
counts[idx[:-6], idx[6:]] = 0
counts[idx[6:], idx[:-6]] = 0
counts[idx[:-7], idx[7:]] = 0
counts[idx[7:], idx[:-7]] = 0

c, l = utils.downsample_resolution(counts, np.array([counts.shape[0]]), 6)

c = mem.cache(iced.normalization.ICE_normalization)(c)
e = mem.cache(get_expected)(c, l)

lengths = fastio.load_lengths("data/ay2013/trophozoites_10000_raw.bed")
counts = fastio.load_counts(
    "data/ay2013/trophozoites_10000_raw.matrix",
    lengths=lengths)

counts = counts.toarray()
counts = counts + counts.T
counts[np.diag_indices_from(counts)] = 0
idx = np.arange(len(counts))

counts[idx[:-1], idx[1:]] = 0
counts[idx[1:], idx[-1:]] = 0

torm = counts.sum(axis=0) == 0
counts[torm] = np.nan
counts[:, torm] = np.nan

counts, lengths = utils.downsample_resolution(counts, lengths, 3)

counts = iced.filter.filter_low_counts(counts, lengths=lengths)
counts = mem.cache(iced.normalization.ICE_normalization)(counts)

torm = counts.sum(axis=0) == 0
counts[torm] = np.nan
counts[:, torm] = np.nan

mapping = mem.cache(get_mapping)(counts, lengths)

counts = counts[lengths.cumsum()[5]:lengths.cumsum()[6],
                lengths.cumsum()[5]:lengths.cumsum()[6]]


lengths = np.array([lengths[6]])
exp = mem.cache(get_expected)(counts, lengths, mapping=mapping)

fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
m = ax.matshow(c / e, cmap="RdBu_r", vmax=2, vmin=0, origin="bottom")
cb = fig.colorbar(m)
cb.set_label("Observed/Expected", fontsize="medium")

ax.set_xticks([])
ax.set_xlabel("Chr 7", fontweight="bold")
ax.set_yticks([])
ax.set_ylabel("Chr 7", fontweight="bold")
ax.set_title("VE", fontweight="bold")
ax.text(-0.05, 1.02, "A", transform=ax.transAxes, fontweight="bold",
        fontsize="large")
fig.savefig("figures/plasmo_ve.pdf")

fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
m = ax.matshow(counts / exp, cmap="RdBu_r", vmax=2, vmin=0, origin="bottom")
cb = fig.colorbar(m)
cb.set_label("Observed/Expected", fontsize="medium")
ax.set_xticks([])
ax.set_xlabel("Chr 7", fontweight="bold")

ax.set_yticks([])
ax.set_ylabel("Chr 7", fontweight="bold")
ax.set_title("HiC", fontweight="bold")
ax.text(-0.05, 1.02, "B", transform=ax.transAxes, fontweight="bold",
        fontsize="large")
fig.savefig("figures/plasmo_hic.pdf")
