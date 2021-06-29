############################################################
##
## DRAGEN Open Source Software
## Copyright (c) 2019-2020 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
############################################################
##
## To configure the build, see config.mk
##
############################################################

include config.mk

all: $(programs:%=$(BUILD)/%)

.PHONY: clean
clean:
	$(RMDIR) $(DRAGEN_OS_BUILD_DIR_BASE)

.PHONY: help
clean:
help: $(DRAGEN_OS_ROOT_DIR)/README.md
	cat $<

############################################################
##
## Rules and includes for the actual build as needed.
## empty MAKECMDGOALS defaults to "all". Inclusion must happen if any goal is not in "clean help"
##
############################################################
ifneq ($(filter-out clean help, $(or $(MAKECMDGOALS), all)),)

# Dependencies are initially generated with ".Td" extension to avoid issues if compiling fails afterwards
# and a POSTCOMPILE operation is needed to rename the file with the final ".d" extension
.PRECIOUS: %.d
DEPFLAGS = -MT $@ -MMD -MP -MF $(@:%.o=%.Td)
POSTCOMPILE ?= mv -f $(@:%.o=%.Td) $(@:%.o=%.d)
%.d: ;

# use a .sentinel file as a proxy to directories to avoid time stamp galore
.PRECIOUS: %/.sentinel
%/.sentinel:
	@mkdir -p $* && touch $@

include $(wildcard $(BUILD)/testRunner.d)

# side effects:
#  - builds 'libraries' variable required for linking programs, integration and system tests
#  - builds and executes unit and integration tests for each librarry
ssw_lib_dirs_aux:=$(SSW_LIBS)
include $(foreach lib_dir, $(SSW_LIBS), $(DRAGEN_OS_MAKE_DIR)/ssw_lib.mk)
dragen_stub_lib_dirs_aux:=$(DRAGEN_STUB_LIBS)
include $(foreach lib_dir, $(DRAGEN_STUB_LIBS), $(DRAGEN_OS_MAKE_DIR)/dragen_stub_lib.mk)
dragen_lib_dirs_aux:=$(DRAGEN_LIBS)
include $(foreach lib_dir, $(DRAGEN_LIBS), $(DRAGEN_OS_MAKE_DIR)/dragen_lib.mk)
lib_dirs_aux:=$(DRAGEN_OS_LIBS)
include $(foreach lib_dir, $(DRAGEN_OS_LIBS), $(DRAGEN_OS_MAKE_DIR)/lib.mk)

programs_aux:=$(programs)
include $(foreach program, $(programs), $(DRAGEN_OS_MAKE_DIR)/program.mk)

# programs for system tests
include $(DRAGEN_OS_MAKE_DIR)/tests.mk

include $(DRAGEN_OS_MAKE_DIR)/install.mk
endif

############################################################
##
## Tracing make variables
## Only add these targets if the goal is to print as it adds
## spurious targets for all non-print goals specified on the command line
##
############################################################
ifneq (,$(filter print-%, $(MAKECMDGOALS)))
print-%: ; @$(error $* is $($*) (from $(origin $*)))
$(filter-out print-%, $(MAKECMDGOALS)): $(filter print-%, $(MAKECMDGOALS))
endif
