from __future__ import print_function
import numpy as np
from glob import glob
from iced import io
from iced import filter
from iced import normalization
import matplotlib.pyplot as plt

from utils import get_mapping


filenames = glob("data/ay2013/*10000_raw.matrix") + \
            glob("data/lemieux2013/25kb/*.matrix")
filenames.sort()


fig, ax = plt.subplots(figsize=(12, 12))

for filename in filenames:
    lengths = io.load_lengths(filename.replace(".matrix", ".bed"))
    counts = io.load_counts(filename, lengths=lengths)

    if "25kb" in filename:
        resolution = 25000
    elif "20000" in filename:
        resolution = 20000
    else:
        resolution = 10000


    counts = counts.toarray()
    counts = counts.T + counts

    # Just making sure there is no interaction counted in teh diag
    counts[np.diag_indices_from(counts)] = 0
    counts = filter.filter_low_counts(counts, percentage=0.03, sparsity=False)
    counts = normalization.ICE_normalization(counts)

    print("1. Compute count vs genomic distance relationship")
    mapping = get_mapping(counts, lengths, verbose=True, smoothed=False)
    ax.plot(mapping[0, 2:] * resolution, mapping[1, 2:])
    ax.axhline(mapping[1, 0])


ax.set_yscale("log")
ax.set_xscale("log")
