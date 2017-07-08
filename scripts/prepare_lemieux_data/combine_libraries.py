from __future__ import print_function
from iced import io
from glob import glob


libraries = ["MobI", "HindIII"]
filenames = glob("data/lemieux/*/*_raw.matrix")
filenames = [filename for filename in filenames
             if not filename.endswith("_combined_raw.matrix")]
filenames.sort()
basenames = set(
    [filename.strip("_raw.matrix").rstrip("MobI").rstrip("HindIII")
     for filename in filenames])
lengths = io.load_lengths("data/ay2013/rings_10000_raw.bed")


for basename in basenames:
    print("Combining", basename)
    counts = None
    for filename in filenames:
        if filename.startswith(basename):
            if counts is None:
                counts = io.load_counts(filename, lengths=lengths)
            else:
                counts += io.load_counts(filename, lengths=lengths)
    io.write_counts(
        basename + "_combined_raw.matrix", counts)
