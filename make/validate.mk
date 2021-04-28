SHELL:=/bin/bash -o pipefail

SAMTOOLS:=samtools

DRAGENOS_DIR:=$(shell pwd)
DRAGEN:=$(DRAGENOS_DIR)/build/release/dragen-os
COMPARE:=$(DRAGENOS_DIR)/build/release/test/compare

FASTQALL:=../FC1-NA12878/original/FC1-NA12878-01_S1.1M.fastq
FASTQR1:=../FC1-NA12878/original/FC1-NA12878-01_S1.R1_1M.fastq
REFERENCE:=../dos-ref
OUTPUT:=../dos
SAM_DIR:=/tmp
ID:=
GOOD_ID:=

STATIC_PREFIX:=mapq60Positionsc-tlsdragen-fastq
BASELINE_PREFIX:=../drgn.flr400/mapq60Positionsc-tlsdragen

PE_STATS:=--Aligner.pe-stat-mean-insert 209.814 --Aligner.pe-stat-stddev-insert 80.7387 \
--Aligner.pe-stat-quartiles-insert '154 198 260' --Aligner.pe-stat-mean-read-len 142.043
STATIC_PARAMS:= --RGID rgid --RGSM rgsm --preserve-map-align-order yes --mmap-reference true --ref-load-hash-bin false
PARAMS:=$(STATIC_PARAMS) $(PE_STATS) -r $(REFERENCE)

default: all

prefix=$(subst $(OUTPUT)/,,$(subst .sam,,$@))
interleaved_aux=$(patsubst $(STATIC_PREFIX)-se-%,no,$(prefix))
interleaved=$(patsubst $(STATIC_PREFIX)-pe-%,yes,$(interleaved_aux))
swall_aux=$(patsubst %-swall-$(ID),1,$(prefix))
swall=$(patsubst %-$(ID),0,$(swall_aux))
fastq_aux=$(patsubst $(STATIC_PREFIX)-se-%,$(FASTQR1),$(prefix))
fastq=$(patsubst $(STATIC_PREFIX)-pe-%,$(FASTQALL),$(fastq_aux))

SAFEPIPETARGET=	"$@.tmp" && mv $@.tmp $@

$(OUTPUT)/$(STATIC_PREFIX)-pe-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-pe-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-se-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-se-$(ID).sam :
	/usr/bin/time -v $(DRAGEN) $(PARAMS) \
	--output-file-prefix $(prefix) -1 $(fastq) --interleaved $(interleaved) --Aligner.sw-all=$(swall) \
	> $(SAFEPIPETARGET)

PE_SAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-pe-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-pe-$(ID).sam

SE_SAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-se-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-se-$(ID).sam

SAM_FILES:=$(PE_SAM_FILES) $(SE_SAM_FILES)

DIFF_FILES:=$(subst .sam,.diff,$(SAM_FILES))
COMP_FILES:=$(subst .sam,.comp,$(SAM_FILES))


$(BASELINE_PREFIX)-%.sam : $(BASELINE_PREFIX)-%.bam
	$(SAMTOOLS) view -h $< > $(SAFEPIPETARGET)

COMPARE_DIR=$(OUTPUT)/compare-$(ID)-$*

$(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp: $(BASELINE_PREFIX)-%.sam $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).sam 
	mkdir -p $(COMPARE_DIR) && (echo -ne "$*\t" && $(COMPARE) $^ -o $(COMPARE_DIR) -k yes --position-mismatches yes |grep -v IMPLEMENTED | grep mismatch |sed 's/=/ /') > $(SAFEPIPETARGET)

$(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).diff: $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).sam
	(echo -ne "$*\t" && (diff <(samtools view $(BASELINE_PREFIX)-$*.bam |sed 's/\tXQ:i:[0-9]*//' |cut -f1,2,3,4,5,6,7,8,9,12,13,14,15 |sed 's/\tSD:f:.*//' ) \
	   <(samtools view $< |sed 's/\tXQ:i:[0-9]*//' |cut -f1,2,3,4,5,6,7,8,9,12,13,14,15); [ $$? -eq 1 ] ) | wc -l) > $(SAFEPIPETARGET)


pe : $(OUTPUT)/$(STATIC_PREFIX)-pe-$(ID).comp
	for f in $^; do head -1 $$f; done

pe-swall : $(OUTPUT)/$(STATIC_PREFIX)-pe-swall-$(ID).comp
	for f in $^; do head -1 $$f; done

se : $(OUTPUT)/$(STATIC_PREFIX)-se-$(ID).comp
	for f in $^; do head -1 $$f; done

se-swall : $(OUTPUT)/$(STATIC_PREFIX)-se-swall-$(ID).comp
	for f in $^; do head -1 $$f; done


$(OUTPUT)/$(STATIC_PREFIX)-compare-%-$(ID).txt : $(OUTPUT)/$(STATIC_PREFIX)-%-$(GOOD_ID).comp $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
	paste -d ' '  $^ |sed 's/$*\t//g' | cut -f 1,2,4 -d' ' > $(SAFEPIPETARGET)

compare-pe : $(OUTPUT)/$(STATIC_PREFIX)-compare-pe-$(ID).txt
	cat $<

compare-pe-swall : $(OUTPUT)/$(STATIC_PREFIX)-compare-pe-swall-$(ID).txt
	cat $<

compare-se : $(OUTPUT)/$(STATIC_PREFIX)-compare-se-$(ID).txt
	cat $<

compare-se-swall : $(OUTPUT)/$(STATIC_PREFIX)-compare-se-swall-$(ID).txt
	cat $<


all: $(COMP_FILES)
	for f in $^; do head -1 $$f; done

clean:
	-rm $(SAM_FILES)
	-rm $(DIFF_FILES)
	-rm $(COMP_FILES)

clean-diffs:
	-rm $(DIFF_FILES)
	-rm $(COMP_FILES)

 