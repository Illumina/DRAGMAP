
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
lib_objects := $(lib_sources:$(SSW_SRC_DIR)/$(lib_dir)/%.cpp=$(BUILD)/ssw_$(lib_dir)/%.o)
lib_objects += $(lib_c_sources:$(SSW_SRC_DIR)/$(lib_dir)/%.c=$(BUILD)/ssw_$(lib_dir)/%.o)

include $(wildcard $(lib_objects:%.o=%.d))

###
# Compile the actual library components and generate the library
#
#
###


$(BUILD)/ssw_$(lib_dir).a: lib_objects:=$(lib_objects)
$(BUILD)/ssw_$(lib_dir).a: $(lib_objects)
	$(AR) crfs $@ $(lib_objects)

#$(BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp $(BUILD)/$(lib_dir)/%.d $(BUILD)/$(lib_dir)/.sentinel
#	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
#	$(POSTCOMPILE)

# Note: the dependency on $(libraries) is to force the order of compilation to be the same as the order of declaration of the libraries
$(BUILD)/ssw_$(lib_dir)/%.o: lib_dir:=$(lib_dir)
$(BUILD)/ssw_$(lib_dir)/%.o: $(SSW_SRC_DIR)/$(lib_dir)/%.cpp $(BUILD)/ssw_$(lib_dir)/%.d $(BUILD)/ssw_$(lib_dir)/.sentinel $(libraries)
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

$(BUILD)/ssw_$(lib_dir)/%.o: $(SSW_SRC_DIR)/$(lib_dir)/%.c $(BUILD)/ssw_$(lib_dir)/%.d $(BUILD)/ssw_$(lib_dir)/.sentinel $(libraries)
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

#$(BUILD)/$(lib_dir)/%.d: ;


###
# must be built in reverse order for linking
###
libraries := $(BUILD)/ssw_$(lib_dir).a $(libraries)
