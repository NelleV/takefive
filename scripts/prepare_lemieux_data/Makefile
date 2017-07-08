DATA = $(wildcard data/lemieux_/*/*.gz)
COUNTS = $(patsubst data/lemieux_/%.gz, data/lemieux/%_raw.matrix, $(DATA))

all: $(COUNTS)

data/lemieux/%_raw.matrix: data/lemieux_/%.gz convert_to_hicpro.py
	python convert_to_hicpro.py $< -o $@
