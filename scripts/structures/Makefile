DATASETS = $(wildcard data/*/*.matrix)
STRUCTURES = $(patsubst data/%.matrix, results/%_UNB0_01_structure.txt, $(DATASETS))
IMG = $(patsubst results/%_structure.txt, images/%_distances.png, $(STRUCTURES))


all: $(STRUCTURES) $(IMG)

results/bunnik2017/Pv_%_UNB0_01_structure.txt: data/bunnik2017/Pv_%.matrix infer_structures_nb.py
	python infer_structures_nb.py $< --lengths data/$*.bed --is-pv

results/%_UNB0_01_structure.txt: data/%.matrix infer_structures_nb.py
	python infer_structures_nb.py $< --lengths data/$*.bed

images/%_UNB0_01_distances.png: results/%_UNB0_01_structure.txt plot_distances.py
	python plot_distances.py $< --lengths data/$*.bed
