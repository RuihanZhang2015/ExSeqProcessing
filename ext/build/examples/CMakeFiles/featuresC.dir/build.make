# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mp/nas1/fixstars/karl/SIFT3D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mp/nas1/fixstars/karl/SIFT3D/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/featuresC.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/featuresC.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/featuresC.dir/flags.make

examples/CMakeFiles/featuresC.dir/featuresC.c.o: examples/CMakeFiles/featuresC.dir/flags.make
examples/CMakeFiles/featuresC.dir/featuresC.c.o: ../examples/featuresC.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/CMakeFiles/featuresC.dir/featuresC.c.o"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/examples && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/featuresC.dir/featuresC.c.o   -c /mp/nas1/fixstars/karl/SIFT3D/examples/featuresC.c

examples/CMakeFiles/featuresC.dir/featuresC.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/featuresC.dir/featuresC.c.i"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/examples && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mp/nas1/fixstars/karl/SIFT3D/examples/featuresC.c > CMakeFiles/featuresC.dir/featuresC.c.i

examples/CMakeFiles/featuresC.dir/featuresC.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/featuresC.dir/featuresC.c.s"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/examples && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mp/nas1/fixstars/karl/SIFT3D/examples/featuresC.c -o CMakeFiles/featuresC.dir/featuresC.c.s

examples/CMakeFiles/featuresC.dir/featuresC.c.o.requires:

.PHONY : examples/CMakeFiles/featuresC.dir/featuresC.c.o.requires

examples/CMakeFiles/featuresC.dir/featuresC.c.o.provides: examples/CMakeFiles/featuresC.dir/featuresC.c.o.requires
	$(MAKE) -f examples/CMakeFiles/featuresC.dir/build.make examples/CMakeFiles/featuresC.dir/featuresC.c.o.provides.build
.PHONY : examples/CMakeFiles/featuresC.dir/featuresC.c.o.provides

examples/CMakeFiles/featuresC.dir/featuresC.c.o.provides.build: examples/CMakeFiles/featuresC.dir/featuresC.c.o


# Object files for target featuresC
featuresC_OBJECTS = \
"CMakeFiles/featuresC.dir/featuresC.c.o"

# External object files for target featuresC
featuresC_EXTERNAL_OBJECTS =

examples/featuresC: examples/CMakeFiles/featuresC.dir/featuresC.c.o
examples/featuresC: examples/CMakeFiles/featuresC.dir/build.make
examples/featuresC: lib/libsift3D.so
examples/featuresC: lib/libimutil.so
examples/featuresC: examples/CMakeFiles/featuresC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable featuresC"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/featuresC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/featuresC.dir/build: examples/featuresC

.PHONY : examples/CMakeFiles/featuresC.dir/build

examples/CMakeFiles/featuresC.dir/requires: examples/CMakeFiles/featuresC.dir/featuresC.c.o.requires

.PHONY : examples/CMakeFiles/featuresC.dir/requires

examples/CMakeFiles/featuresC.dir/clean:
	cd /mp/nas1/fixstars/karl/SIFT3D/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/featuresC.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/featuresC.dir/clean

examples/CMakeFiles/featuresC.dir/depend:
	cd /mp/nas1/fixstars/karl/SIFT3D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mp/nas1/fixstars/karl/SIFT3D /mp/nas1/fixstars/karl/SIFT3D/examples /mp/nas1/fixstars/karl/SIFT3D/build /mp/nas1/fixstars/karl/SIFT3D/build/examples /mp/nas1/fixstars/karl/SIFT3D/build/examples/CMakeFiles/featuresC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/featuresC.dir/depend
