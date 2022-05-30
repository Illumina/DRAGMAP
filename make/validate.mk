SHELL:=/bin/bash -o pipefail

SAMTOOLS:=samtools

USE_BAM ?= 0
DRAGENOS_DIR:=$(shell pwd)
DRAGEN:=$(DRAGENOS_DIR)/build/release/dragen-os
COMPARE:=$(DRAGENOS_DIR)/build/release/compare

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
STATIC_PARAMS:= --RGID rgid --RGSM rgsm --ref-load-hash-bin false --preserve-map-align-order true
PARAMS:=$(STATIC_PARAMS) $(PE_STATS) -r $(REFERENCE) $(USER_PARAMS)

default: all

prefix=$(subst .sam,,$(subst .bam,,$(subst tmp.,,$(subst $(OUTPUT)/,,$@))))
interleaved_aux=$(patsubst $(STATIC_PREFIX)-se-%,,$(prefix))
inter=$(patsubst $(STATIC_PREFIX)-pe-%,--interleaved,$(interleaved_aux))
swall_aux=$(patsubst %-swall-$(ID),1,$(prefix))
swall=$(patsubst %-$(ID),0,$(swall_aux))
fastq_aux=$(patsubst $(STATIC_PREFIX)-se-%,$(FASTQR1),$(prefix))
fastq=$(patsubst $(STATIC_PREFIX)-pe-%,$(FASTQALL),$(fastq_aux))

SAFEPIPETARGET=	"$@.tmp" && mv $@.tmp $@

PE_SAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-pe-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-pe-$(ID).sam

SE_SAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-se-swall-$(ID).sam \
$(OUTPUT)/$(STATIC_PREFIX)-se-$(ID).sam

SAM_FILES:=$(PE_SAM_FILES) $(SE_SAM_FILES)

ifneq (1,$(USE_BAM))
$(SAM_FILES) :
	/usr/bin/time -v $(DRAGEN) $(PARAMS) \
	--output-file-prefix $(prefix) -1 $(fastq) $(inter) --Aligner.sw-all=$(swall) \
	> $(SAFEPIPETARGET)
else

PE_BAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-pe-swall-$(ID).bam \
$(OUTPUT)/$(STATIC_PREFIX)-pe-$(ID).bam

SE_BAM_FILES:= $(OUTPUT)/$(STATIC_PREFIX)-se-swall-$(ID).bam \
$(OUTPUT)/$(STATIC_PREFIX)-se-$(ID).bam

BAM_FILES:=$(PE_BAM_FILES) $(SE_BAM_FILES)

$(BAM_FILES) :
	mkdir -p $(OUTPUT) && \
	LC_ALL=C /usr/bin/time -v $(DRAGEN) $(PARAMS) \
	--output-directory $(OUTPUT) --output-file-prefix tmp.$(prefix) -1 $(fastq) $(inter) --Aligner.sw-all=$(swall) && \
	rename tmp.$(prefix) $(prefix) $(OUTPUT)/tmp.$(prefix)*
endif

COMP_FILES:=$(subst .sam,.comp,$(SAM_FILES))

#.PRECIOUS: $(BASELINE_PREFIX)-%.sam
#$(BASELINE_PREFIX)-%.sam : $(BASELINE_PREFIX)-%.bam
#	$(SAMTOOLS) view -h $< > $(SAFEPIPETARGET)

.PRECIOUS: %.sam
%.sam : %.bam
	$(SAMTOOLS) view -h $< > $(SAFEPIPETARGET)

COMPARE_DIR=$(OUTPUT)/compare-$(ID)-$*

.PRECIOUS: $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
$(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp: $(BASELINE_PREFIX)-%.sam $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).sam 
	mkdir -p $(COMPARE_DIR) && (echo -ne "$*\t" && $(COMPARE) $^ -o $(COMPARE_DIR) -k yes --position-mismatches yes |grep -v IMPLEMENTED |grep -v _SD_ | grep mismatch |sed 's/=/ /') > $(SAFEPIPETARGET)

% : $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
	head -1 $<

# .comp files. Do not require GOOD_ID run. Just produces the diffs against the baseline bam.
comp: pe pe-swall se se-swall

compare-% : $(BASELINE_PREFIX)-%-$(GOOD_ID).comp $(OUTPUT)/$(STATIC_PREFIX)-%-$(ID).comp
	paste -d ' '  $(foreach f,$^,<(sort $f)) |sed 's/$*\t//g' | cut -f 1,2,4 -d' '

# Use these to see if the new build is an improvement over the GOOD_ID one
compare : compare-se-swall compare-pe-swall compare-se compare-pe

sam: $(SAM_FILES)

baseline: $(SAM_FILES)
	rename -- '-$(ID)' '' $(SAM_FILES)

all: $(COMP_FILES)
	for f in $^; do head -1 $$f; done

clean:
	-rm $(SAM_FILES)
	-rm $(BAM_FILES)
	-rm $(COMP_FILES)

clean-diffs:
	-rm $(COMP_FILES)

