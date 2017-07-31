import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from minorswing import fastio
import iced
from joblib import Memory


mem = Memory(".joblib")


fig, axes = plt.subplots(ncols=2, figsize=(8, 4), tight_layout=True)

# Load troph data
lengths = fastio.load_lengths(
    "data/lemieux2013/25kb/A4_pilot_combined_raw.bed")

counts = fastio.load_counts(
    "data/lemieux2013/25kb/A4_pilot_combined_raw.matrix",
    lengths=lengths)
c = fastio.load_counts(
    "data/lemieux2013/25kb/A4+_combined_raw.matrix",
    lengths=lengths)
counts += c
c = fastio.load_counts(
    "data/lemieux2013/25kb/A44_combined_raw.matrix",
    lengths=lengths)
counts += c

centromeres = np.loadtxt("files/pf.cent").mean(axis=1)
centromeres /= 25000

counts[np.diag_indices_from(counts)] = 0
counts = iced.filter.filter_low_counts(counts, sparsity=False)
counts = mem.cache(iced.normalization.ICE_normalization)(counts)

counts = counts.toarray()
counts = counts.T + counts

# Extract chromosome 7
begin, end = lengths.cumsum()[5], lengths.cumsum()[6]
counts = counts[begin:end, begin:end]
length = lengths[6]

to_rm = counts.sum(axis=0) == 0
counts[np.diag_indices_from(counts)] = np.nanmax(counts)
counts[to_rm] = np.nan
counts[:, to_rm] = np.nan

ax = axes[0]
ax.text(-0.05, 1.02, "A", transform=ax.transAxes, fontweight="bold",
fontsize="large")

ax.matshow(
    counts,
    norm=colors.SymLogNorm(3),
    cmap="YlGnBu",
    extent=(0, length, 0, length),
    origin="bottom")

ax.axhline(centromeres[6], linestyle="-", color="#000000", alpha=0.5)
ax.axvline(centromeres[6], linestyle="-", color="#000000", alpha=0.5)

ax.set_xlabel("Chr 7", fontweight="bold", fontsize="small")
ax.set_ylabel("Chr 7", fontweight="bold", fontsize="small")
ax.set_title("Lemieux-A44/A4/A+", fontweight="bold", fontsize="medium")

ax.set_xticks([])
ax.set_yticks([])

##########################################
# Load troph data
lengths = fastio.load_lengths(
    "data/ay2013/trophozoites_10000_raw.bed")
counts = fastio.load_counts(
    "data/ay2013/trophozoites_10000_raw.matrix",
    lengths=lengths)

centromeres = np.loadtxt("files/pf.cent").mean(axis=1)
centromeres /= 10000

counts = iced.filter.filter_low_counts(counts, sparsity=False)
counts = mem.cache(iced.normalization.ICE_normalization)(counts)

counts = counts.toarray()
counts = counts.T + counts
to_rm = counts.sum(axis=0) == 0
counts[to_rm] = np.nan
counts[:, to_rm] = np.nan

# Extract chromosome 7
begin, end = lengths.cumsum()[5], lengths.cumsum()[6]
counts = counts[begin:end, begin:end]
length = lengths[6]

ax = axes[1]
ax.text(-0.05, 1.02, "B", transform=ax.transAxes, fontweight="bold",
fontsize="large")
ax.matshow(
    counts,
    norm=colors.SymLogNorm(0.05),
    cmap="YlGnBu",
    extent=(0, length, 0, length),
    vmax=14000,
    origin="bottom")
ax.axhline(centromeres[6], linestyle="-", color="#000000", alpha=0.5)
ax.axvline(centromeres[6], linestyle="-", color="#000000", alpha=0.5)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("Chr 7", fontweight="bold", fontsize="small")
ax.set_ylabel("Chr 7", fontweight="bold", fontsize="small")
ax.set_title("Ay-trophozoite", fontweight="bold", fontsize="medium")

try:
    os.makedirs("figures")
except OSError:
    pass


fig.savefig("figures/counts_maps.png")
