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

lengths = io.load_lengths("data/ay2013/rings_10000_raw.bed")
counts = np.zeros([lengths.sum(), lengths.sum()])

begin_chr = np.concatenate([[0], lengths.cumsum()])

with gzip.open(filename, "rb") as f:
    for i, line in enumerate(f):
        chrom1, loc1, chrom2, loc2, c = line.split("\t")
        loc1, loc2, c = int(loc1), int(loc2), int(c)

        chrom1 = int(chrom1.strip("chr"))
        chrom2 = int(chrom2.strip("chr"))

        loc1 /= 10000
        loc2 /= 10000

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
