
ifeq (,$(ssw_lib_dirs_aux))
$(error No library directories specified)
endif

lib_dir := $(word 1, $(ssw_lib_dirs_aux))
ifeq (1,$(words $(ssw_lib_dirs_aux)))
ssw_lib_dirs_aux:=
else
ssw_lib_dirs_aux:=$(wordlist 2, $(words $(ssw_lib_dirs_aux)), $(ssw_lib_dirs_aux))
endif

lib_sources := $(wildcard $(SSW_SRC_DIR)/$(lib_dir)/*.cpp)
lib_c_sources := $(wildcard $(SSW_SRC_DIR)/$(lib_dir)/*.c)
lib_objects := $(lib_sources:$(SSW_SRC_DIR)/$(lib_dir)/%.cpp=$(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.o)
lib_objects += $(lib_c_sources:$(SSW_SRC_DIR)/$(lib_dir)/%.c=$(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.o)

include $(wildcard $(lib_objects:%.o=%.d))

###
# Compile the actual library components and generate the library
#
#
###


$(DRAGEN_OS_BUILD)/libdragmap-$(lib_dir).a: lib_objects:=$(lib_objects)
$(DRAGEN_OS_BUILD)/libdragmap-$(lib_dir).a: $(lib_objects)
	$(SILENT) $(AR) crfs $@ $(lib_objects)

#$(DRAGEN_OS_BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp $(DRAGEN_OS_BUILD)/$(lib_dir)/%.d $(DRAGEN_OS_BUILD)/$(lib_dir)/.sentinel
#	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
#	$(POSTCOMPILE)

# Note: the dependency on $(libraries) is to force the order of compilation to be the same as the order of declaration of the libraries
$(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.o: lib_dir:=$(lib_dir)
$(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.o: $(SSW_SRC_DIR)/$(lib_dir)/%.cpp $(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.d $(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/.sentinel
	$(SILENT) $(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $< && $(POSTCOMPILE)

$(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.o: $(SSW_SRC_DIR)/$(lib_dir)/%.c $(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/%.d $(DRAGEN_OS_BUILD)/ssw_$(lib_dir)/.sentinel
	$(SILENT) $(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $< && $(POSTCOMPILE)

#$(DRAGEN_OS_BUILD)/$(lib_dir)/%.d: ;


###
# must be built in reverse order for linking
###
libraries := $(DRAGEN_OS_BUILD)/libdragmap-$(lib_dir).a $(libraries)

###
# targets for external tools
###
library_targets := $(lib_dir)-lib $(library_targets)
.PHONY: $(lib_dir)-lib
$(lib_dir)-lib : $(DRAGEN_OS_BUILD)/libdragmap-$(lib_dir).a


