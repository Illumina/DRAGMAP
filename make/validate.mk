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

#USER_PARAMS:= --Aligner.vectorized-sw yes
PE_STATS:=--Aligner.pe-stat-mean-insert 209.814 --Aligner.pe-stat-stddev-insert 80.7387 \
--Aligner.pe-stat-quartiles-insert '154 198 260' --Aligner.pe-stat-mean-read-len 142.043
STATIC_PARAMS:= --RGID rgid --RGSM rgsm --mmap-reference true --ref-load-hash-bin false --preserve-map-align-order yes 
PARAMS:=$(STATIC_PARAMS) $(PE_STATS) -r $(REFERENCE) $(USER_PARAMS)

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

COMP_FILES:=$(subst .sam,.comp,$(SAM_FILES))

.PRECIOUS: $(BASELINE_PREFIX)-%.sam
$(BASELINE_PREFIX)-%.sam : $(BASELINE_PREFIX)-%.bam
	$(SAMTOOLS) view -h $< > $(SAFEPIPETARGET)

COMPARE_DIR=$(OUTPUT)/compare-$(ID)-$*

.PRECIOUS: $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
$(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp: $(BASELINE_PREFIX)-%.sam $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).sam 
	mkdir -p $(COMPARE_DIR) && (echo -ne "$*\t" && $(COMPARE) $^ -o $(COMPARE_DIR) -k yes --position-mismatches yes |grep -v IMPLEMENTED | grep mismatch |sed 's/=/ /') > $(SAFEPIPETARGET)

% : $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
	head -1 $<

compare-% : $(OUTPUT)/$(STATIC_PREFIX)-%-$(GOOD_ID).comp $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
	paste -d ' '  $^ |sed 's/$*\t//g' | cut -f 1,2,4 -d' '

compare : compare-se-swall compare-pe-swall compare-se compare-pe

all: $(COMP_FILES)
	for f in $^; do head -1 $$f; done

clean:
	-rm $(SAM_FILES)
	-rm $(COMP_FILES)

clean-diffs:
	-rm $(COMP_FILES)

