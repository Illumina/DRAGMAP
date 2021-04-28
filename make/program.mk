ifeq (,$(programs_aux))
$(error No programs specified)
endif

program := $(word 1, $(programs_aux))

ifeq (1, $(words $(programs_aux)))
programs_aux:=
else
programs_aux:=$(wordlist 2, $(words $(programs_aux)), $(programs_aux))
endif

$(BUILD)/$(program).o: $(DRAGEN_OS_SRC_DIR)/$(program).cpp $(BUILD)/$(program).d $(BUILD)/.sentinel
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

$(BUILD)/$(program): $(BUILD)/$(program).o $(libraries)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(libraries) $(LDFLAGS)

#$(BUILD)/$(program).d: ;

include $(wildcard $(BUILD)/$(program).d)

