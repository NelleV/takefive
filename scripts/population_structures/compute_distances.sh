python compute_downsampled_distances.py \
  structures/ay2013/rings_10000_raw_PO \
  --lengths data/ay2013/rings_10000_raw.bed

python compute_downsampled_distances.py \
  structures/ay2013/schizonts_10000_raw_PO \
  --lengths data/ay2013/rings_10000_raw.bed

python compute_downsampled_distances.py \
  structures/ay2013/trophozoites_10000_raw_PO \
  --lengths data/ay2013/rings_10000_raw.bed

python compute_downsampled_distances.py \
  structures/lemieux2013/25kb/B15C2_combined_raw_PO \
  --factor 4 \
  --lengths data/lemieux2013/25kb/DCJ_On_combined_raw.bed
