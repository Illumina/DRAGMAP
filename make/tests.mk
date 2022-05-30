# building all available system tests

TEST_BUILD_DIR=$(DRAGEN_OS_BUILD)/test

ifeq (1,$(HAS_GTEST))
system_tests:=$(patsubst $(DRAGEN_OS_TEST_DIR)/%.cpp,%,$(wildcard $(DRAGEN_OS_TEST_DIR)/*Gtest.cpp))
system_tests:=$(filter-out $(system_skipped), $(system_tests))

define SYSTEM_TEST

test_programs: $(system_tests:%=$(TEST_BUILD_DIR)/%)
all: test_programs

system_tool := $(1)

$(TEST_BUILD_DIR)/$(1).o: $(DRAGEN_OS_TEST_DIR)/$(1).cpp $(TEST_BUILD_DIR)/$(1).d $(TEST_BUILD_DIR)/.sentinel
	$(SILENT_SE) $$(CXX) $$(DEPFLAGS) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -c -o $$@ $$< && $$(POSTCOMPILE)

$(TEST_BUILD_DIR)/$(1): $(TEST_BUILD_DIR)/$(1).o $(libraries)
	$(SILENT_SE) $$(CXX) $$(CPPFLAGS) $$(GTEST_CPPFLAGS) $$(CXXFLAGS) -o $$@ $$< $$(libraries)  $(GTEST_LDFLAGS) -lgtest_main -lgtest $$(LDFLAGS)

#$(DRAGEN_OS_BUILD)/system/$(system_tool).d: ;
include $(wildcard $(TEST_BUILD_DIR)/$(1).d)

endef # define SYSTEM_TEST

$(foreach t,$(system_tests),$(eval $(call SYSTEM_TEST,$(t))))
endif

# building tools that are independent of gtest

system_tools:=$(patsubst $(DRAGEN_OS_TEST_DIR)/%.cpp,%,$(wildcard $(DRAGEN_OS_TEST_DIR)/*.cpp))
system_tools:=$(filter-out %Gtest, $(system_tools))

tools_programs: $(system_tools:%=$(TEST_BUILD_DIR)/%)
all: tools_programs

define SYSTEM_TOOL

system_tool := $(1)

$(TEST_BUILD_DIR)/$(1).o: $(DRAGEN_OS_TEST_DIR)/$(1).cpp $(TEST_BUILD_DIR)/$(1).d $(TEST_BUILD_DIR)/.sentinel
	$(SILENT_SE) $$(CXX) $$(DEPFLAGS) $$(CPPFLAGS) $$(CXXFLAGS) -c -o $$@ $$< && $$(POSTCOMPILE)

$(TEST_BUILD_DIR)/$(1): $(TEST_BUILD_DIR)/$(1).o $(libraries)
	$(SILENT_SE) $$(CXX) $$(CPPFLAGS) $$(CXXFLAGS) -o $$@ $$< $$(libraries) $$(LDFLAGS)

#$(DRAGEN_OS_BUILD)/system/$(system_tool).d: ;
include $(wildcard $(TEST_BUILD_DIR)/$(1).d)

endef # define SYSTEM_TOOL

$(foreach t,$(system_tools),$(eval $(call SYSTEM_TOOL,$(t))))

