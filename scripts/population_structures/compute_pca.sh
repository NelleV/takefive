
for chr in {1..14}; do
  python clustering_across_stages_downsampled_distances.py \
    -c $chr

  python clustering_across_stages_downsampled_distances.py \
    -c $chr -a UNB0
done;

python clustering_across_stages_downsampled_distances.py
python clustering_across_stages_downsampled_distances.py -a UNB0
python clustering_across_stages_downsampled_distances.py -a PO

