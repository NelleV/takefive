import argparse
import numpy as np
from scipy import sparse

from iced import io
from iced import utils

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--outname", "-o")
args = parser.parse_args()

filename = args.filename
lengths = io.load_lengths("data/ay2013/rings_10000_raw.bed")


counts = io.load_counts(filename, lengths=lengths)
counts = counts.toarray()
counts = counts.T + counts

new_counts, new_lengths = utils.downsample_resolution(counts, lengths)

new_counts = sparse.coo_matrix(np.triu(new_counts))


io.write_counts(args.outname, new_counts)
io.write_lengths(args.outname.replace(".matrix", ".bed"), new_lengths)
