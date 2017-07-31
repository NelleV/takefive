import numpy as np
import argparse
from iced import io
import utils
from scipy import spatial
from statistics import compute_mepd, compute_witten_and_noble
import matplotlib.pyplot as plt
from sklearn.externals.joblib import Memory

mem = Memory(".joblib")

fig, axes = plt.subplots(nrows=4) 

filename = "structures/lemieux2013/25kb/B15C2_combined_raw_PO_01_structure.txt"
lengths = io.load_lengths("data/lemieux2013/25kb/B15C2_combined_raw.bed")
#filename = "structures/ay2013/trophozoites_10000_raw_PO_01_structure.txt"
#lengths = io.load_lengths("data/ay2013/trophozoites_10000_raw.bed")
X = np.loadtxt(filename)
centromeres = np.loadtxt("files/pf.cent")
chr_ = np.arange(len(centromeres))

if "25kb" in filename:
    resolution = 25000
elif "20000" in filename:
    resolution = 20000
else:
    resolution = 10000

centromeres = np.round(centromeres.mean(axis=1) / resolution)
centromeres = centromeres.astype(int)

cent_pval_l, cent_pval_r, cent_mepd, cent_null = mem.cache(compute_mepd)(
    X, lengths, centromeres, chr_,
    return_null=True)
ax = axes[0]
ax.hist(cent_null, bins=100, color="#000000")
ax.axvline(cent_mepd, color="#AB0000")

cent_pval_l, cent_pval_r, cent_mepd, cent_null = mem.cache(
    compute_witten_and_noble)(
        X, lengths, centromeres, chr_,
        threshold=0.2,
        return_null=True)
ax = axes[1]
nbins = len(np.unique(cent_null))
ax.hist(cent_null, bins=nbins, color="#000000")
ax.axvline(cent_mepd, color="#AB0000")

# Telomeres:
telomeres = np.concatenate(
    [[0], lengths[:-1], lengths - 1])
chr_ = np.tile(chr_, 2)

tel_pval_l, tel_pval_r, tel_mepd, tel_null = compute_mepd(
    X, lengths, telomeres, chr_,
    return_null=True)

