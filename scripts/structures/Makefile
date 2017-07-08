DATASETS = $(wildcard data/ay2013/*.matrix) \
	   $(wildcard data/lemieux2013/25kb/*.matrix)

STRUCTURES = $(patsubst data/%.matrix, results/%_UNB0_01_structure.txt, $(DATASETS)) \
	     $(patsubst data/%.matrix, results/%_PO_01_structure.txt, $(DATASETS))  \
	     $(patsubst data/%.matrix, results/%_MDS_01_structure.txt, $(DATASETS)) 

IMG = $(patsubst results/%_structure.txt, images/%_distances.png, $(STRUCTURES))


all: $(STRUCTURES) $(IMG)

results/bunnik2017/Pv_%_UNB0_01_structure.txt: data/bunnik2017/Pv_%.matrix infer_structures_nb.py
	python infer_structures_nb.py $< --lengths data/$*.bed --is-pv

results/%_UNB0_01_structure.txt: data/%.matrix infer_structures_nb.py
	python infer_structures_nb.py $< --lengths data/$*.bed

images/%_UNB0_01_distances.png: results/%_UNB0_01_structure.txt plot_distances.py
	python plot_distances.py $< --lengths data/$*.bed


results/%_MDS_01_structure.txt: data/%.matrix infer_structures_mds.py
	python infer_structures_mds.py $< --lengths data/$*.bed

images/%_MDS_01_distances.png: results/%_MDS_01_structure.txt plot_distances.py
	python plot_distances.py $< --lengths data/$*.bed

results/%_PO_01_structure.txt: data/%.matrix infer_structures_po.py
	python infer_structures_po.py $< --lengths data/$*.bed

images/%_PO_01_distances.png: results/%_PO_01_structure.txt plot_distances.py
	python plot_distances.py $< --lengths data/$*.bed
