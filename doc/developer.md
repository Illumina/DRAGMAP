
# Directory structure

## Top level

Only the strict minimum:

* README.md
* LICENSE and COPYRIGHT
* Makefile
* config.mk

## make subdirectory

All the included makefiles for the build system

## doc subdirectory

All aspects of the documentation:

* usage
* developer
* method

## make subdirectory

All the included components for the top level Makefile.

## src

All the C++ source code. The headers are under include and the compilation units under lib.
Th include and src subdirectories have the same top level structure: one subdirectory dor each
library.

Each of the subdirectory under lib should have a 'tests' subditectory for all the unit
tests for the corresponding library. See the "Unit tests" section below for details.

The 'tests' subdirectory only contains the source code for running the tests.

# Coding conventions

As a general rule, a single class ("SomeClass") has a corresponding declaration file ("SomeClass.hh")
in a subdirectory of include and a definition file ("SomeClass.cpp") in the same subdirectory of lib.
See "src/include/options/DragenOsOptions.hh" and "src/lib/options/DragenOsOptions.cpp".
In some cases it makes sense to aggregate multiple classes into a single file.
See "src/include/common/Exceptions.hh" and "src/lib/common/Exceptions.cpp" for examplars.

Style preference: currently there is a mixture of many difference styles but
the preference should be towards camelCase vesus underscore_separated with:
* UpperCamelCase for Classes
* lowerCamelCase for methods, functions, parameter names and local variables
* m_UpperCameLCase for class members

# Build system

## Overview

The build system is a modular (using includes) non-recursive Makefile.
The Makefile initially includes the configuration "config.mk". It then includes "program.mk" for each of
the target programs and "lib.mk" for each of the library - which in turns includes "gtest.mk" for the unit
tests for the corresponding library. Finally, it includes "integration.mk" for each of the integration
tests and "install.mk" just once.

The lists of programs, libraries, compilation units (CU) and integration tests are discovered with glob expressions:

* programs: "*.cpp"
* libraries: "src/lib/*/*.cpp" then using the unique paths. The first wiltcard is used for the library name. The
  second wildcard is used for the list of compilation units for each library.
* unit tests: "src/lib/$(library)/tests/*Fixture.cpp"
* integration tests: "tests/integration/*.cpp"

The build output is either the subdirectory "build/release" or "build/debug", depending on the build type.

Dependencies are auto-built while compiling each CU and included upon availability.

There is some basic analysis of MAKECMDGOALS to avoid all the inclusions - and particularly the generation of
the dependencies - when the list of goals explicitly contains only generis targets like "clean" or "help".

## Programs

Compilation with DEPFLAGS, CPPFLAGS and CXXFLAGS. Linking with CPPFLAGS, CXXFLAGS and LDFLAGS. There is a dummy
empty target for the dependency file for each program.  There is an include of the corresponding dependency file
for each program.

DEPFLAGS are just the usual macro "-MT $@ -MMD -MP -MF $(@:%.o=%.Td)" to generate the dependencies. Note that
the dependency files are generated with the extension ".Td" as a temporary file and then renamed with the ".d"
extension in a POSTCOMPILE instruction - this avoids dependency issues when the compile crashes after the
generation of the dependencies.

CPPFLAGS specify c++11, all warnings, add boost as additional system include path, optimization and debug level.

CXXFLAGS controls the outpu of the compiler depending on the build type. In debug mode, it would only adding the 
sanitize (address, leak and undefined), coverage and profiling information. In release mode, it would be all the
meaningful optimizations.

LDFLAGS specifies the path to the Boost libraries, the list of system libraries and a few additional linker
flags like '--exclude-libraries=ALL' and such.

## Libraries

The list of expected libraries must be explicitly defined in the correct order for linking (built static). As
a sanity check, the makefile verified that the list of expected libraries matches the list of libraries found
with the glob expression for that library. Library CUs are compiles as the programs are with DEPFLAGS, CPPFLAGS
and CXXFLAGS and a POSTCOMPILE to rename the dependency file. There is also an include for the dependency file
for each CU.

Each library include "make/gtest.mk" for the optional unit tests associated to each CU.

## Unit tests

The unit tests are built and run in isolation for each CU. The main program to execute the unit tests is 
the default provided in libgteest_main.

The unit tests are for each CU - they are optional though. For a CU CompilationUnit, three files are expected:
"CompilationUnitMocks.hh", "CompilationUnitFixture.hh" and "CompilationUnitFixture.cpp". The build system first
generate a wrapper "CompilationUnitWrapper.cpp":

    #include "CompilationUnitMocks.hh"
    #include "CompilationUnit.cpp"

Both "CompilationUnitFixture.cpp" and "CompilationUnitWrapper.cpp" are compiled and then linked to "testRunner.o"
to create "testCompilationUnit" which is then executed to create "passedCompilationUnit" which is then added to 
the list of dependencies for the library.

# Testing

## Integration tests

Each integration test is compiled and linked as the programs are. 

TODO: add support for the automation of the execution of these tests.

## Acceptance tests

Not supported

