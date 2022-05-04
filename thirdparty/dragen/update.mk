FIND:=find
GREP:=grep
CP:=cp

ifeq (,$(DRAGEN_SRC))
$(error Please set DRAGEN_SRC to the location of dragen repository root folder)
endif

DRAGMAP_FILES:=$(shell $(FIND) . -type f |$(GREP) -v $(lastword $(MAKEFILE_LIST)) | $(GREP) -v mapping_stats)

all: $(DRAGMAP_FILES)

.SECONDEXPANSION:
$(DRAGMAP_FILES) : $(DRAGEN_SRC)/$$@
	$(CP) $< $@.tmp && \
	sed -i '/^\#ifndef OPEN_SOURCE/,/^\#endif  \/\/ OPEN_SOURCE/{d}' $@.tmp && \
	mv $@.tmp $@ \

