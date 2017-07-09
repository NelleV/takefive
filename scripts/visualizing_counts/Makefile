COUNTS = $(wildcard data/ay2013/*_raw.matrix) \
	 $(wildcard data/lemieux2013/*/*combined_raw.matrix) \

OE = $(patsubst data/%.matrix, results/%_OE.matrix, $(COUNTS)) \
	
IMG = $(patsubst data/%.matrix, images/%.png, $(COUNTS)) \
      $(patsubst results/%.matrix, images/%.png, $(OE))

all: $(IMG) $(OE)

images/%_raw.png: data/%_raw.matrix plot_counts.py
	python plot_counts.py $< --lengths data/$*_raw.bed

images/%_OE.png: results/%_OE.matrix plot_OE.py
	python plot_OE.py $< --lengths data/$*.bed -o $@

results/%_raw_OE.matrix: data/%_raw.matrix compute_observed_over_expected.py
	python compute_observed_over_expected.py $< -o $@ --bed-file data/$*_raw.bed --normalize
