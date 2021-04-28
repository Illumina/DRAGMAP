
DRAGEN_OS_FLAGS:=-include linux/limits.h

ifeq (,$(dragen_lib_dirs_aux))
$(error No library directories specified)
endif

lib_dir := $(word 1, $(dragen_lib_dirs_aux))
ifeq (1,$(words $(dragen_lib_dirs_aux)))
dragen_lib_dirs_aux:=
else
dragen_lib_dirs_aux:=$(wordlist 2, $(words $(dragen_lib_dirs_aux)), $(dragen_lib_dirs_aux))
endif

lib_sources := $(wildcard $(DRAGEN_SRC_DIR)/$(lib_dir)/*.cpp)
lib_c_sources := $(wildcard $(DRAGEN_SRC_DIR)/$(lib_dir)/*.c)
lib_objects := $(lib_sources:$(DRAGEN_SRC_DIR)/$(lib_dir)/%.cpp=$(BUILD)/drgn_$(lib_dir)/%.o)
lib_objects += $(lib_c_sources:$(DRAGEN_SRC_DIR)/$(lib_dir)/%.c=$(BUILD)/drgn_$(lib_dir)/%.o)

include $(wildcard $(lib_objects:%.o=%.d))

###
# Compile the actual library components and generate the library
#
#
###
$(BUILD)/drgn_$(lib_dir).a: lib_objects:=$(lib_objects)
$(BUILD)/drgn_$(lib_dir).a: $(lib_objects)
	$(AR) crfs $@ $(lib_objects)

#$(BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp $(BUILD)/$(lib_dir)/%.d $(BUILD)/$(lib_dir)/.sentinel
#	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
#	$(POSTCOMPILE)

# Note: the dependency on $(libraries) is to force the order of compilation to be the same as the order of declaration of the libraries
$(BUILD)/drgn_$(lib_dir)/%.o: lib_dir:=$(lib_dir)
$(BUILD)/drgn_$(lib_dir)/%.o: $(DRAGEN_SRC_DIR)/$(lib_dir)/%.cpp $(BUILD)/drgn_$(lib_dir)/%.d $(BUILD)/drgn_$(lib_dir)/.sentinel $(libraries)
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(DRAGEN_OS_FLAGS) -c -o $@ $<
	$(POSTCOMPILE)

$(BUILD)/drgn_$(lib_dir)/%.o: $(DRAGEN_SRC_DIR)/$(lib_dir)/%.c $(BUILD)/drgn_$(lib_dir)/%.d $(BUILD)/drgn_$(lib_dir)/.sentinel $(libraries)
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) $(DRAGEN_OS_FLAGS) -c -o $@ $<
	$(POSTCOMPILE)

#$(BUILD)/$(lib_dir)/%.d: ;


###
# must be built in reverse order for linking
###
libraries := $(BUILD)/drgn_$(lib_dir).a $(libraries)
