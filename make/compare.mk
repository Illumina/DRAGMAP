SHELL:=/bin/bash

DRAGENOS_DIR:=$(shell pwd)
BASELINE_SAM:=../drgn.flr400/mapq60Positionsc-tlsdragen.sam
FASTQ1:=../FC1-NA12878/original/FC1-NA12878-01_S1.R1_1M.fastq
FASTQ2:=../FC1-NA12878/original/FC1-NA12878-01_S1.R2_1M.fastq
#REFERENCE=/ssd/dragenos/reference_genomes/Hsapiens/hg38-noalt-with-decoy/DRAGEN/8/
REFERENCE:=../dos-ref
OUTPUT:=/tmp
SAM_DIR:=/tmp
ID:=$(shell git rev-parse HEAD | cut -b1-8)

default: all

SAM:=$(OUTPUT)/new-$(ID)-FC1-NA12878-01_S1.R1_R2_1M.sam

$(SAM):
	echo $(ID) ; date ; /usr/bin/time -v $(DRAGENOS_DIR)/build/release/dragen-os \
	-1 "$(FASTQ1)" -2 "$(FASTQ2)" -r "$(REFERENCE)"  \
	--preserve-map-align-order 1 --Aligner.pe-orientation 0 \
	--Aligner.pe-stat-mean-insert 210 --Aligner.pe-stat-mean-read-len 148 --Aligner.pe-stat-quartiles-insert "154 199 262" \
	--Aligner.pe-stat-stddev-insert 82 --Aligner.rescue-ceil-factor 5 \
	--enable-sampling 0 --ref-load-hash-bin 1 --mmap-reference 1 > "$@.tmp" 2> "$@.out" && mv $@.tmp $@

COMPARE_TXT:=/tmp/compare-$(ID).txt
$(COMPARE_TXT): $(SAM)
	$(DRAGENOS_DIR)/build/release/test/compare $(BASELINE_SAM) "$<" -k 1 > $@

all: /tmp/compare-$(ID).txt
	paste <(cat concordance.txt) <(cut -d= -f2 /tmp/compare-$(ID).txt) | awk -e '{print $$1 "|" $$2 "|"}'

clean:
	-rm $(SAM)
	-rm $(COMPARE_TXT)

