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
include imutil/CMakeFiles/meximutil.dir/depend.make

# Include the progress variables for this target.
include imutil/CMakeFiles/meximutil.dir/progress.make

# Include the compile flags for this target's objects.
include imutil/CMakeFiles/meximutil.dir/flags.make

imutil/CMakeFiles/meximutil.dir/imutil.c.o: imutil/CMakeFiles/meximutil.dir/flags.make
imutil/CMakeFiles/meximutil.dir/imutil.c.o: ../imutil/imutil.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object imutil/CMakeFiles/meximutil.dir/imutil.c.o"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/meximutil.dir/imutil.c.o   -c /mp/nas1/fixstars/karl/SIFT3D/imutil/imutil.c

imutil/CMakeFiles/meximutil.dir/imutil.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/meximutil.dir/imutil.c.i"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mp/nas1/fixstars/karl/SIFT3D/imutil/imutil.c > CMakeFiles/meximutil.dir/imutil.c.i

imutil/CMakeFiles/meximutil.dir/imutil.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/meximutil.dir/imutil.c.s"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mp/nas1/fixstars/karl/SIFT3D/imutil/imutil.c -o CMakeFiles/meximutil.dir/imutil.c.s

imutil/CMakeFiles/meximutil.dir/imutil.c.o.requires:

.PHONY : imutil/CMakeFiles/meximutil.dir/imutil.c.o.requires

imutil/CMakeFiles/meximutil.dir/imutil.c.o.provides: imutil/CMakeFiles/meximutil.dir/imutil.c.o.requires
	$(MAKE) -f imutil/CMakeFiles/meximutil.dir/build.make imutil/CMakeFiles/meximutil.dir/imutil.c.o.provides.build
.PHONY : imutil/CMakeFiles/meximutil.dir/imutil.c.o.provides

imutil/CMakeFiles/meximutil.dir/imutil.c.o.provides.build: imutil/CMakeFiles/meximutil.dir/imutil.c.o


imutil/CMakeFiles/meximutil.dir/nifti.c.o: imutil/CMakeFiles/meximutil.dir/flags.make
imutil/CMakeFiles/meximutil.dir/nifti.c.o: ../imutil/nifti.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object imutil/CMakeFiles/meximutil.dir/nifti.c.o"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/meximutil.dir/nifti.c.o   -c /mp/nas1/fixstars/karl/SIFT3D/imutil/nifti.c

imutil/CMakeFiles/meximutil.dir/nifti.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/meximutil.dir/nifti.c.i"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mp/nas1/fixstars/karl/SIFT3D/imutil/nifti.c > CMakeFiles/meximutil.dir/nifti.c.i

imutil/CMakeFiles/meximutil.dir/nifti.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/meximutil.dir/nifti.c.s"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mp/nas1/fixstars/karl/SIFT3D/imutil/nifti.c -o CMakeFiles/meximutil.dir/nifti.c.s

imutil/CMakeFiles/meximutil.dir/nifti.c.o.requires:

.PHONY : imutil/CMakeFiles/meximutil.dir/nifti.c.o.requires

imutil/CMakeFiles/meximutil.dir/nifti.c.o.provides: imutil/CMakeFiles/meximutil.dir/nifti.c.o.requires
	$(MAKE) -f imutil/CMakeFiles/meximutil.dir/build.make imutil/CMakeFiles/meximutil.dir/nifti.c.o.provides.build
.PHONY : imutil/CMakeFiles/meximutil.dir/nifti.c.o.provides

imutil/CMakeFiles/meximutil.dir/nifti.c.o.provides.build: imutil/CMakeFiles/meximutil.dir/nifti.c.o


imutil/CMakeFiles/meximutil.dir/dicom.cpp.o: imutil/CMakeFiles/meximutil.dir/flags.make
imutil/CMakeFiles/meximutil.dir/dicom.cpp.o: ../imutil/dicom.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object imutil/CMakeFiles/meximutil.dir/dicom.cpp.o"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/meximutil.dir/dicom.cpp.o -c /mp/nas1/fixstars/karl/SIFT3D/imutil/dicom.cpp

imutil/CMakeFiles/meximutil.dir/dicom.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/meximutil.dir/dicom.cpp.i"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mp/nas1/fixstars/karl/SIFT3D/imutil/dicom.cpp > CMakeFiles/meximutil.dir/dicom.cpp.i

imutil/CMakeFiles/meximutil.dir/dicom.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/meximutil.dir/dicom.cpp.s"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mp/nas1/fixstars/karl/SIFT3D/imutil/dicom.cpp -o CMakeFiles/meximutil.dir/dicom.cpp.s

imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.requires:

.PHONY : imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.requires

imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.provides: imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.requires
	$(MAKE) -f imutil/CMakeFiles/meximutil.dir/build.make imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.provides.build
.PHONY : imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.provides

imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.provides.build: imutil/CMakeFiles/meximutil.dir/dicom.cpp.o


# Object files for target meximutil
meximutil_OBJECTS = \
"CMakeFiles/meximutil.dir/imutil.c.o" \
"CMakeFiles/meximutil.dir/nifti.c.o" \
"CMakeFiles/meximutil.dir/dicom.cpp.o"

# External object files for target meximutil
meximutil_EXTERNAL_OBJECTS =

lib/wrappers/matlab/libmeximutil.so: imutil/CMakeFiles/meximutil.dir/imutil.c.o
lib/wrappers/matlab/libmeximutil.so: imutil/CMakeFiles/meximutil.dir/nifti.c.o
lib/wrappers/matlab/libmeximutil.so: imutil/CMakeFiles/meximutil.dir/dicom.cpp.o
lib/wrappers/matlab/libmeximutil.so: imutil/CMakeFiles/meximutil.dir/build.make
lib/wrappers/matlab/libmeximutil.so: /usr/local/MATLAB/R2018a/bin/glnxa64/libmex.so
lib/wrappers/matlab/libmeximutil.so: /usr/local/MATLAB/R2018a/bin/glnxa64/libmwlapack.so
lib/wrappers/matlab/libmeximutil.so: /usr/local/MATLAB/R2018a/bin/glnxa64/libmwblas.so
lib/wrappers/matlab/libmeximutil.so: /usr/lib64/libz.so
lib/wrappers/matlab/libmeximutil.so: /usr/lib64/libm.so
lib/wrappers/matlab/libmeximutil.so: imutil/CMakeFiles/meximutil.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mp/nas1/fixstars/karl/SIFT3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared library ../lib/wrappers/matlab/libmeximutil.so"
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/meximutil.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
imutil/CMakeFiles/meximutil.dir/build: lib/wrappers/matlab/libmeximutil.so

.PHONY : imutil/CMakeFiles/meximutil.dir/build

imutil/CMakeFiles/meximutil.dir/requires: imutil/CMakeFiles/meximutil.dir/imutil.c.o.requires
imutil/CMakeFiles/meximutil.dir/requires: imutil/CMakeFiles/meximutil.dir/nifti.c.o.requires
imutil/CMakeFiles/meximutil.dir/requires: imutil/CMakeFiles/meximutil.dir/dicom.cpp.o.requires

.PHONY : imutil/CMakeFiles/meximutil.dir/requires

imutil/CMakeFiles/meximutil.dir/clean:
	cd /mp/nas1/fixstars/karl/SIFT3D/build/imutil && $(CMAKE_COMMAND) -P CMakeFiles/meximutil.dir/cmake_clean.cmake
.PHONY : imutil/CMakeFiles/meximutil.dir/clean

imutil/CMakeFiles/meximutil.dir/depend:
	cd /mp/nas1/fixstars/karl/SIFT3D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mp/nas1/fixstars/karl/SIFT3D /mp/nas1/fixstars/karl/SIFT3D/imutil /mp/nas1/fixstars/karl/SIFT3D/build /mp/nas1/fixstars/karl/SIFT3D/build/imutil /mp/nas1/fixstars/karl/SIFT3D/build/imutil/CMakeFiles/meximutil.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : imutil/CMakeFiles/meximutil.dir/depend
