COUNTS = $(wildcard data/ay2013/*10000*.matrix)
SCOUNTS = $(patsubst %_10000_raw.matrix, %_20000_raw.matrix, $(COUNTS))


all: $(SCOUNTS)

%_20000_raw.matrix: %_10000_raw.matrix subsample_20kb.py
	python subsample_20kb.py $< -o $@

