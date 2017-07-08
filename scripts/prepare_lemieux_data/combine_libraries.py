from __future__ import print_function
from iced import io
from glob import glob


libraries = ["MobI", "HindIII"]
filenames = glob("data/lemieux/*/*_raw.matrix")
filenames = [filename for filename in filenames
             if not filename.endswith("_combined_raw.matrix")]
filenames.sort()
basenames = list(set(
    [filename.split("-")[0]
     for filename in filenames]))
basenames.sort()

for basename in basenames:
    print("Combining", basename)
    counts = None
    for filename in filenames:
        if filename.startswith(basename):

            lengths = io.load_lengths(filename.replace("_raw.matrix", ".bed"))
            if counts is None:
                counts = io.load_counts(filename, lengths=lengths)
            else:
                counts += io.load_counts(filename, lengths=lengths)
    io.write_counts(
        basename + "_combined_raw.matrix", counts)
    io.write_lengths(
        basename + "_combined_raw.bed", lengths)
