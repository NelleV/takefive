from __future__ import print_function
import numpy as np
from glob import glob
from iced import io
from iced import utils

filenames = glob("data/ay2013/*10000_raw.matrix") + \
            glob("data/lemieux2013/25kb/*.matrix")
filenames.sort()

for filename in filenames:
    lengths = io.load_lengths(filename.replace(".matrix", ".bed"))
    counts = io.load_counts(filename, lengths=lengths)

    counts = counts.toarray()
    counts = counts.T + counts

    mask = utils.get_intra_mask(lengths)

    # Just making sure there is no interaction counted in teh diag
    counts[np.diag_indices_from(counts)] = 0
    print(filename)
    print("Total number of counts", counts.sum())
    print("%% of intra", counts[mask].sum()/counts.sum() * 100)
    print("%% of inter", counts[np.invert(mask)].sum()/counts.sum() * 100)
    print()
