# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /local/home/xiaoyang1/anaconda3/envs/lgn_data/bin/cmake

# The command to remove a file.
RM = /local/home/xiaoyang1/anaconda3/envs/lgn_data/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /local/home/xiaoyang1/HEPData4ML-master/util/root

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build

# Utility rule file for move_header_vectorcalcs.

# Include the progress variables for this target.
include vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/progress.make

vectorcalcs/CMakeFiles/move_header_vectorcalcs: include/vectorcalcs/VectorCalcs.h


include/vectorcalcs/VectorCalcs.h: ../inc/vectorcalcs/VectorCalcs.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Copying header /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/inc/vectorcalcs/VectorCalcs.h to /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/include"
	cd /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/vectorcalcs && /local/home/xiaoyang1/anaconda3/envs/lgn_data/bin/cmake -E copy /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/inc/vectorcalcs/VectorCalcs.h /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/include/vectorcalcs/VectorCalcs.h

move_header_vectorcalcs: vectorcalcs/CMakeFiles/move_header_vectorcalcs
move_header_vectorcalcs: include/vectorcalcs/VectorCalcs.h
move_header_vectorcalcs: vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/build.make

.PHONY : move_header_vectorcalcs

# Rule to build all files generated by this target.
vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/build: move_header_vectorcalcs

.PHONY : vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/build

vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/clean:
	cd /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/vectorcalcs && $(CMAKE_COMMAND) -P CMakeFiles/move_header_vectorcalcs.dir/cmake_clean.cmake
.PHONY : vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/clean

vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/depend:
	cd /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /local/home/xiaoyang1/HEPData4ML-master/util/root /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/vectorcalcs /local/home/xiaoyang1/HEPData4ML-master/util/root/vectorcalcs/build/vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : vectorcalcs/CMakeFiles/move_header_vectorcalcs.dir/depend

