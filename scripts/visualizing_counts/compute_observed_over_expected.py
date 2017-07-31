import os
import argparse

import io_ as io
from iced import filter
from iced import normalization
from utils import get_mapping, get_expected


parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--normalize", "-n", action="store_true", default=False)
parser.add_argument("--bed-file", "-b")
parser.add_argument("--outfile", "-o")
args = parser.parse_args()

lengths, base = io.load_lengths(args.bed_file, return_base=True)
counts = io.load_counts(args.filename, lengths=lengths, base=base)

if args.normalize:
    counts = filter.filter_low_counts(counts, percentage=0.03, sparsity=False)
    counts = normalization.ICE_normalization(counts)

print("1. Compute count vs genomic distance relationship")
mapping = get_mapping(counts, lengths, verbose=True)

print("2. Estimating expected...")
c_expected = get_expected(counts, lengths, mapping=mapping)

print("3. Estimating observed over expected...")
counts.data /= c_expected


if args.outfile is not None:
    try:
        os.makedirs(os.path.dirname(args.outfile))
    except OSError:
        pass
    io.write_counts(args.outfile, counts)
