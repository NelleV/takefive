COUNTS = $(wildcard data/ay2013/*_raw.matrix) \
	 $(wildcard data/lemieux2013/*/*combined_raw.matrix) \

IMG = $(patsubst data/%.matrix, images/%.png, $(COUNTS))

all: $(IMG)

images/%.png: data/%.matrix plot_counts.py
	python plot_counts.py $< --lengths data/$*.bed

