
ifeq (,$(lib_dirs_aux))
$(error No library directories specified)
endif

lib_dir := $(word 1, $(lib_dirs_aux))
ifeq (1,$(words $(lib_dirs_aux)))
lib_dirs_aux:=
else
lib_dirs_aux:=$(wordlist 2, $(words $(lib_dirs_aux)), $(lib_dirs_aux))
endif

lib_sources := $(wildcard $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/*.cpp)
lib_c_sources := $(wildcard $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/*.c)
lib_objects := $(lib_sources:$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp=$(BUILD)/$(lib_dir)/%.o)
lib_objects += $(lib_c_sources:$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.c=$(BUILD)/$(lib_dir)/%.o)

include $(wildcard $(lib_objects:%.o=%.d))

###
# Build and run the unit tests before actually building the library
# Strict googletest unit tests (testing of a single compilation unit)
#
# In practice, each compilation unit (CU) in the library should have a
# unique test fixture (TF). However, as unit testing is not enforced,
# we have to account for untested CUs.
# TODO: enforce unit testing for all CUs
#
# Method: 
# -------
# Each Compilation Unit (CU) is tested independently. Header-only components
# can be tested as if there was a completely empty corresponding CU.
# the relevant files in the $(unit_test_src_dir) are:
# - $(CU)Gtest.cpp: the gtest tests for the CU
# - $(CU)Mocks.hpp: the mocks to use for the CU - if any.
# Each CU is compiled with the corresponding mocks - if any - to enable
# deep testing of system libraries (e.g. testing vector memory allocation
# failures or boost filesystem errors).
#
# Note: these tests should NOT have any dependency on the "system", particularly
# no IOs (no file system, no network, rtc.)
#
# Note: each CU is wrapped in a corresponding $(CU)Wrapper.cpp that also 
# include the optional mocks before the source code for the CU. This could have
# been achieved with the g++ command line option "-include $(CU)Mocks.hpp".
# Thhe choice for a frapper file is driven by two factors: (1) support for 
# header-only CUs and (2) clarity - otherwise, the unit test $(CU).o could
# be misleading for a casual observer.
###

unit_test_src_dir:=$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/tests/unit
unit_test_build_dir:=$(BUILD)/$(lib_dir)/tests/unit

unit_test_sources:=$(wildcard $(unit_test_src_dir)/*Gtest.cpp)
unit_tests:=$(unit_test_sources:$(unit_test_src_dir)/%Gtest.cpp=%)

define UNIT_GTEST
include $(wildcard $(unit_test_build_dir)/$(1)Wrapper.d $(unit_test_build_dir)/$(1)Gtest.d)
.PRECIOUS: $(unit_test_build_dir)/$(1)Wrapper.cpp $(unit_test_build_dir)/$(1)Wrapper.d 
$(unit_test_build_dir)/$(1)Wrapper.cpp: $(unit_test_build_dir)/.sentinel
	$(ECHO) > $$@ ; \
	[[ -e $(unit_test_src_dir)/$(1)Mocks.hpp ]] && $(ECHO) \#include \"$(1)Mocks.hpp\" >> $$@ ; \
	[[ -e $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/$(1).cpp ]] && $(ECHO) \#include \"$(1).cpp\" >> $$@ ; \
	$(ECHO) >> $$@ # important for case without .cpp

$(unit_test_build_dir)/$(1)Wrapper.o: $(unit_test_build_dir)/$(1)Wrapper.cpp $(unit_test_build_dir)/$(1)Wrapper.d
	$$(CXX) -I$(unit_test_src_dir) -I$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir) $$(DEPFLAGS) $$(CPPFLAGS) $$(CXXFLAGS) -c -o $$@ $$<
	$$(POSTCOMPILE)

$(unit_test_build_dir)/$(1)Gtest.o: $(unit_test_src_dir)/$(1)Gtest.cpp $(unit_test_build_dir)/$(1)Gtest.d $(unit_test_build_dir)/.sentinel
	$(CXX) -I$(unit_test_src_dir) -I$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir) $$(DEPFLAGS) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -c -o $$@ $$<
	$$(POSTCOMPILE)

$(unit_test_build_dir)/$(1)Gtest: $(unit_test_build_dir)/$(1)Gtest.o $(unit_test_build_dir)/$(1)Wrapper.o
	$$(CXX) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -o $$@ $$^ $(GTEST_LDFLAGS) $$(LDFLAGS)

# For developers convenience, shows the output of failed tests
$(unit_test_build_dir)/$(1)Gtest.passed: $(unit_test_build_dir)/$(1)Gtest
	$(SILENT){ \
	TEST_COMMAND="$$< > $$<.failed && $$(MV) $$<.failed $$@" ; \
	$(ECHO) "$$$${TEST_COMMAND}" ; \
	$(EVAL) "$$$${TEST_COMMAND}"  ; \
	} || { \
	rc=$$$$? ; $(ECHO) ; $(ECHO) '***' ERROR '***' Failed unit test $$< "($(unit_test_src_dir)/$(1)Gtest.cpp)" ; \
	$(ECHO) ; $(ECHO) $$<.failed: ; $(ECHO) ; $(CAT) $$<.failed ; \
	exit $$$$rc ; \
	}

$$(BUILD)/$(lib_dir).a: $(unit_test_build_dir)/$(1)Gtest.passed

endef # UNIT_GTEST

ifeq (1,${HAS_GTEST})
$(foreach t, $(unit_tests), $(eval $(call UNIT_GTEST,$(t))))
endif

###
# Compile the actual library components and generate the library
#
#
###
$(BUILD)/$(lib_dir).a: lib_objects:=$(lib_objects)
$(BUILD)/$(lib_dir).a: $(lib_objects)
	$(AR) crfs $@ $(lib_objects)

#$(BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp $(BUILD)/$(lib_dir)/%.d $(BUILD)/$(lib_dir)/.sentinel
#	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
#	$(POSTCOMPILE)

# Note: the dependency on $(libraries) is to force the order of compilation to be the same as the order of declaration of the libraries
$(BUILD)/$(lib_dir)/%.o: lib_dir:=$(lib_dir)
$(BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.cpp $(BUILD)/$(lib_dir)/%.d $(BUILD)/$(lib_dir)/.sentinel $(libraries)
	$(CXX) $(DEPFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

$(BUILD)/$(lib_dir)/%.o: $(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/%.c $(BUILD)/$(lib_dir)/%.d $(BUILD)/$(lib_dir)/.sentinel $(libraries)
	$(CC) $(DEPFLAGS) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<
	$(POSTCOMPILE)

#$(BUILD)/$(lib_dir)/%.d: ;

###
# build and run the integration tests for the library
# googletest integration tests (interaction between compilation units)
#
# In practice, each compilation unit (CU) in the library should have
# test fixtures (TF). However, as integration testing is not enforced,
# we have to account for untested CUs.
# TODO: enforce integration testing for all CUs
# TODO: add support for mocks - would require to substitute in lib_objects the
# CU with mocks (either using a wrapper as with the unit tests above or using
# the g++ commandline option "-include $(CU)Mocks.hpp") instead of the original
# CU.
#
# Method: 
# -------
# Each Compilation Unit (CU) is tested in the context of all the compilation units
# built so far (all previous libraries and all CUs from the current library). In
# principle these shouldn't need mocks. In practice, it is likely that some of the
# components won't be designed in a way that enable appropriate level of testing.
# The relevant files in the $(integration_test_src_dir) are:
# - $(CU)Gtest.cpp: the gtest tests for the CU
# - $(CU)Mocks.hpp: the mocks to use for the CU - if any.
#
# Note: these tests should NOT have any dependency on the "system", particularly
# no IOs (no file system, no network, rtc.)
###

integration_test_src_dir:=$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir)/tests/integration
integration_test_build_dir:=$(BUILD)/$(lib_dir)/tests/integration

integration_test_sources:=$(wildcard $(integration_test_src_dir)/*Gtest.cpp)
integration_tests:=$(integration_test_sources:$(integration_test_src_dir)/%Gtest.cpp=%)

define INTEGRATION_GTEST
include $(wildcard $(integration_test_build_dir)/$(1)Gtest.d)

$(integration_test_build_dir)/$(1)Gtest.o: $(integration_test_src_dir)/$(1)Gtest.cpp $(integration_test_build_dir)/$(1)Gtest.d $(integration_test_build_dir)/.sentinel
	$(CXX) -I$(integration_test_src_dir) -I$(DRAGEN_OS_SRC_DIR)/lib/$(lib_dir) $$(DEPFLAGS) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -c -o $$@ $$<
	$$(POSTCOMPILE)

$(integration_test_build_dir)/$(1)Gtest: libraries:= $(libraries)
$(integration_test_build_dir)/$(1)Gtest: lib_objects:=$(lib_objects)
$(integration_test_build_dir)/$(1)Gtest: $(integration_test_build_dir)/$(1)Gtest.o $(lib_objects) $(libraries)
	$$(CXX) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -o $$@ $$^ $(GTEST_LDFLAGS) $$(GTEST_LDFLAGS) $$(LDFLAGS)

# For developers convenience, shows the output of failed tests
$(integration_test_build_dir)/$(1)Gtest.passed: $(integration_test_build_dir)/$(1)Gtest
	$(SILENT){ \
	TEST_COMMAND="$$< > $$<.failed && $$(MV) $$<.failed $$@" ; \
	$(ECHO) "$$$${TEST_COMMAND}" ; \
	$(EVAL) "$$$${TEST_COMMAND}"  ; \
	} || { \
	rc=$$$$? ; $(ECHO) ; $(ECHO) '***' ERROR '***' Failed integration test $$< "($(integration_test_src_dir)/$(1)Gtest.cpp)" ; \
	$(ECHO) ; $(ECHO) $$<.failed: ; $(ECHO) ; $(CAT) $$<.failed ; \
	exit $$$$rc ; \
	}

$$(BUILD)/$(lib_dir).a: $(integration_test_build_dir)/$(1)Gtest.passed

endef # INTEGRATION_GTEST

ifeq (1,${HAS_GTEST})
$(foreach t, $(integration_tests), $(eval $(call INTEGRATION_GTEST,$(t))))
endif


###
# must be built in reverse order for linking
###
libraries := $(BUILD)/$(lib_dir).a $(libraries)

