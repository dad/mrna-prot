PY = $$PYTHONEXE

DATA_DIR = ../data
DATA = ../data
DROPBOX = c:/Dropbox\ \(Drummond\ Lab\)
LIBSRC = ../../lib/src
RDATA = ../../data
RDATA_SCER = $(RDATA)/scerevisiae
GENOME = $(DROPBOX)/bio/genomes
SCER_GEN = $(GENOME)/S288C_reference_genome_Current_Release


features:
	$(PY) make-yeast-seq-features.py $(SCER_GEN)/orf_coding_all_R64-2-1_20150113.fasta $(SCER_GEN)/orf_trans_all_R64-2-1_20150113.fasta $(DATA)/bg/SGD_features.tab $(DATA)/bg/ygob-fungal-pillars-v7.txt --out $(DATA)/bg/scerevisiae-coding-features.txt