# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = "/Applications/CMake 2.8-12.app/Contents/bin/cmake"

# The command to remove a file.
RM = "/Applications/CMake 2.8-12.app/Contents/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = "/Applications/CMake 2.8-12.app/Contents/bin/ccmake"

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /usr/local/include/eigen3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /usr/local/include/eigen3/build_dir

# Include any dependencies generated for this target.
include doc/examples/CMakeFiles/class_Block.dir/depend.make

# Include the progress variables for this target.
include doc/examples/CMakeFiles/class_Block.dir/progress.make

# Include the compile flags for this target's objects.
include doc/examples/CMakeFiles/class_Block.dir/flags.make

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o: doc/examples/CMakeFiles/class_Block.dir/flags.make
doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o: ../doc/examples/class_Block.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /usr/local/include/eigen3/build_dir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o"
	cd /usr/local/include/eigen3/build_dir/doc/examples && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/class_Block.dir/class_Block.cpp.o -c /usr/local/include/eigen3/doc/examples/class_Block.cpp

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/class_Block.dir/class_Block.cpp.i"
	cd /usr/local/include/eigen3/build_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /usr/local/include/eigen3/doc/examples/class_Block.cpp > CMakeFiles/class_Block.dir/class_Block.cpp.i

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/class_Block.dir/class_Block.cpp.s"
	cd /usr/local/include/eigen3/build_dir/doc/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /usr/local/include/eigen3/doc/examples/class_Block.cpp -o CMakeFiles/class_Block.dir/class_Block.cpp.s

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.requires:
.PHONY : doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.requires

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.provides: doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.requires
	$(MAKE) -f doc/examples/CMakeFiles/class_Block.dir/build.make doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.provides.build
.PHONY : doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.provides

doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.provides.build: doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o

# Object files for target class_Block
class_Block_OBJECTS = \
"CMakeFiles/class_Block.dir/class_Block.cpp.o"

# External object files for target class_Block
class_Block_EXTERNAL_OBJECTS =

doc/examples/class_Block: doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o
doc/examples/class_Block: doc/examples/CMakeFiles/class_Block.dir/build.make
doc/examples/class_Block: doc/examples/CMakeFiles/class_Block.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable class_Block"
	cd /usr/local/include/eigen3/build_dir/doc/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/class_Block.dir/link.txt --verbose=$(VERBOSE)
	cd /usr/local/include/eigen3/build_dir/doc/examples && ./class_Block >/usr/local/include/eigen3/build_dir/doc/examples/class_Block.out

# Rule to build all files generated by this target.
doc/examples/CMakeFiles/class_Block.dir/build: doc/examples/class_Block
.PHONY : doc/examples/CMakeFiles/class_Block.dir/build

doc/examples/CMakeFiles/class_Block.dir/requires: doc/examples/CMakeFiles/class_Block.dir/class_Block.cpp.o.requires
.PHONY : doc/examples/CMakeFiles/class_Block.dir/requires

doc/examples/CMakeFiles/class_Block.dir/clean:
	cd /usr/local/include/eigen3/build_dir/doc/examples && $(CMAKE_COMMAND) -P CMakeFiles/class_Block.dir/cmake_clean.cmake
.PHONY : doc/examples/CMakeFiles/class_Block.dir/clean

doc/examples/CMakeFiles/class_Block.dir/depend:
	cd /usr/local/include/eigen3/build_dir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /usr/local/include/eigen3 /usr/local/include/eigen3/doc/examples /usr/local/include/eigen3/build_dir /usr/local/include/eigen3/build_dir/doc/examples /usr/local/include/eigen3/build_dir/doc/examples/CMakeFiles/class_Block.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/examples/CMakeFiles/class_Block.dir/depend

