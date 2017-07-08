import numpy as np
import os
from scipy import sparse
from iced import io
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--outname", "-o")
args = parser.parse_args()

filename = args.filename


lengths = np.loadtxt("files/lengths").astype(float)
if "10kb" in filename:
    resolution = 10000
else:
    resolution = 25000

lengths /= resolution
lengths = np.ceil(lengths).astype(int)


counts = np.zeros([lengths.sum(), lengths.sum()])

begin_chr = np.concatenate([[0], lengths.cumsum()])

with gzip.open(filename, "rb") as f:
    for i, line in enumerate(f):
        chrom1, loc1, chrom2, loc2, c = line.split("\t")
        loc1, loc2, c = int(loc1), int(loc2), int(c)

        chrom1 = int(chrom1.strip("chr"))
        chrom2 = int(chrom2.strip("chr"))

        loc1 /= resolution
        loc2 /= resolution

        loc1 += begin_chr[chrom1 - 1]
        loc2 += begin_chr[chrom2 - 1]

        counts[loc1, loc2] += c

counts = counts.T + counts

if args.outname is not None:
    try:
        os.makedirs(os.path.dirname(args.outname))
    except OSError:
        pass

    io.write_counts(args.outname, sparse.coo_matrix(np.triu(counts)))
    io.write_lengths(args.outname.replace("_raw.matrix", ".bed"), lengths)
