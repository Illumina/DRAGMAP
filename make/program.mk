ifeq (,$(programs_aux))
$(error No programs specified)
endif

program := $(word 1, $(programs_aux))

ifeq (1, $(words $(programs_aux)))
programs_aux:=
else
programs_aux:=$(wordlist 2, $(words $(programs_aux)), $(programs_aux))
endif

$(DRAGEN_OS_BUILD)/$(program).o: $(DRAGEN_OS_SRC_DIR)/$(program).cpp $(DRAGEN_OS_BUILD)/$(program).d $(DRAGEN_OS_BUILD)/.sentinel
	$(SILENT) $(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $< && $(POSTCOMPILE)

$(DRAGEN_OS_BUILD)/$(program): $(DRAGEN_OS_BUILD)/$(program).o $(libraries)
	$(SILENT) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(libraries) $(LDFLAGS)

#$(DRAGEN_OS_BUILD)/$(program).d: ;

include $(wildcard $(DRAGEN_OS_BUILD)/$(program).d)

