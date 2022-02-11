# Required env variables:
# - LBL_PASSWORD
# - LBL_USER

# MAM 20210928 added data/ANI
.PHONY: prepare
prepare: data/ANI mirror_ANIs mirror_unicycler mirror_seqs mirror_genomes

.PHONY: align
align: build/alignments.csv

build data data/ANI data/reference_genomes:
	mkdir -vp $@

# 1. Download ANIm percentage identity from LBL (/data/iarpa/qc/pyani/[SAMPLE]/pyani_out/ANIm_percentage_identity.tab)
# MAM added "mirror_ANIs:" lines below ".PHONY: mirror_ANIs" and ./move_anis.sh
.PHONY: mirror_ANIs
mirror_ANIs:
	sshpass -p ${LBL_PASSWORD} rsync -avR ${LBL_USER}@merlot.lbl.gov:/data/iarpa/qc/pyani/*/pyani_out/ANIm_percentage_identity.tab .
	./src/move_anis.sh

# 2. Download unicycler_assembly genomes from LBL (/data/iarpa/TE/[SAMPLE]/assembly/unicycler_assembly.fasta)
.PHONY: mirror_unicycler
mirror_unicycler:
	sshpass -p ${LBL_PASSWORD} rsync -avRL ${LBL_USER}@merlot.lbl.gov:/data/iarpa/TE/*/assembly/unicycler_assembly.fasta .

# 3. Select best match to unicycler_assembly
build/reference_genomes.csv: src/pick_genomes.py | build
	python3 $< data/ANI $@

# 4. Download best match reference genomes from LBL (/data/iarpa/TE/[SAMPLE]/reference_genome/[GENOME_ID])
# WARNING: this will take quite some time!
# MAM 20210928 sed pattern ^([0-9]+),(.+)$$ can't handle "samples" like 245-original
# deplete them somehow, or modify pattern?
.PHONY: mirror_genomes
mirror_genomes: build/reference_genomes.csv | data/reference_genomes
	$(eval CP=$(shell tail -n +2 build/reference_genomes.csv | sed -E 's|^([0-9]+[^,]*),(.+)$$|/data/iarpa/TE/\1/reference_genome/\2.f*|g' | tr '\r\n' ' '))
	sshpass -p ${LBL_PASSWORD} scp -r -v -T ${LBL_USER}@merlot.lbl.gov:"$(CP)" data/reference_genomes

# 5. Download the target sequences to blast against the genomes
# --recursive --copy-links
# -avR
.PHONY: mirror_seqs
mirror_seqs:
	sshpass -p ${LBL_PASSWORD} rsync --recursive --copy-links --relative ${LBL_USER}@merlot.lbl.gov:/data/iarpa/TE/*/target_sequence/IF*.fasta .

# 6. BLAST! This will also take some time - add '-v' to the command to see progress
build/alignments.csv: src/blast_seqs.py build/reference_genomes.csv
	python3 $^ $@

clean:
	rm -rf build/*
	rm -rf data/*

