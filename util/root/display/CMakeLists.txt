# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.1...3.3)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
  EventDisplay
  VERSION 1.0
  LANGUAGES C CXX)

# Set C++ standard to match ROOT/Delphes
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the path to this repository -- relative to this CMake file.
# NOTE: If this executable is made into a separate git submodule,
#       some of this will likely have to change.
#       Maybe have external installations modify CMAKE_PREFIX_PATH?
set(REPO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

#######################################################
#########
# ROOT
#########

# Find ROOT - as of ~6.20, this sets up some CMake macros too.
find_package(ROOT CONFIG REQUIRED)
#######################################################

#######################################################
#########
# DELPHES
#########

# Find DELPHES.
# Will search CMAKE_PREFIX_PATH, plus the DELPHES_EXT_DIR
# that we define here (to pick up a possible local installation).
file(GLOB DELPHES_SUBDIRS "${REPO_DIR}/external/delphes/*")
foreach(subdir ${DELPHES_SUBDIRS})
    if(IS_DIRECTORY ${subdir})
        set(DELPHES_EXT_DIR ${subdir})
        break()
    endif()
endforeach()

# Check if we found a directory
if(NOT DELPHES_EXT_DIR)
    set(DELPHES_EXT_DIR "${REPO_DIR}/external/delphes/delphes") # TODO: Make this better -- catch the error elegantly
    # message(FATAL_ERROR "No Delphes subdirectory found in ${REPO_DIR}/external/delphes/")
endif()

# TODO: Not sure if this should be handled differently -- I think we want things from
#       classes, modules and display.

# Step 1: Search first for local installation. Using NO_DEFAULT_PATH.
find_path(DELPHES_INCLUDE_DIR
    NAMES classes/DelphesClasses.h
    PATHS
        ${DELPHES_EXT_DIR}
        ${DELPHES_EXT_DIR}/include
    PATH_SUFFIXES delphes Delphes
    NO_DEFAULT_PATH
)

# Step 2: If the above doesn't work, manually set the path.
if(NOT DELPHES_INCLUDE_DIR)
    # Direct path assignment as fallback
    if(EXISTS "${DELPHES_EXT_DIR}/classes/DelphesClasses.h")
        set(DELPHES_INCLUDE_DIR "${DELPHES_EXT_DIR}")
    elseif(EXISTS "${DELPHES_EXT_DIR}/include/classes/DelphesClasses.h")
        set(DELPHES_INCLUDE_DIR "${DELPHES_EXT_DIR}/include")
    else()
        message(FATAL_ERROR "classes/DelphesClasses.h not found in ${DELPHES_EXT_DIR} or ${DELPHES_EXT_DIR}/include")
    endif()
endif()

# Find the "default" Delphes library.
find_library(DELPHES_LIBRARIES
    NAMES Delphes libDelphes
    PATHS
        ${DELPHES_EXT_DIR}/lib
        ${DELPHES_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT DELPHES_LIBRARIES)
    find_library(DELPHES_LIBRARIES
        NAMES Delphes libDelphes # searches CMAKE_PREFIX_PATH
    )
endif()

# Find the Delphes visualization library.
find_library(DELPHES_DISPLAY_LIBRARIES
    NAMES DelphesDisplay libDelphesDisplay
    PATHS
        ${DELPHES_EXT_DIR}/lib
        ${DELPHES_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT DELPHES_DISPLAY_LIBRARIES)
    find_library(DELPHES_DISPLAY_LIBRARIES
        NAMES DelphesDisplay libDelphesDisplay # searches CMAKE_PREFIX_PATH
    )
endif()

# Check if Delphes was found
if(NOT DELPHES_INCLUDE_DIR)
    message(FATAL_ERROR "Delphes headers not found. Please set DELPHES_EXT_DIR environment variable or CMAKE_PREFIX_PATH")
endif()

if(NOT DELPHES_LIBRARIES)
    message(FATAL_ERROR "Delphes library not found. Please set DELPHES_EXT_DIR environment variable or CMAKE_PREFIX_PATH")
endif()

if(NOT DELPHES_DISPLAY_LIBRARIES)
    message(FATAL_ERROR "DelphesDisplay library not found. Please set DELPHES_EXT_DIR environment variable or CMAKE_PREFIX_PATH")
endif()


message(STATUS "Found Delphes headers: ${DELPHES_INCLUDE_DIR}")
message(STATUS "Found Delphes library: ${DELPHES_LIBRARIES}")
message(STATUS "Found DelphesDisplay library: ${DELPHES_DISPLAY_LIBRARIES}")

#######################################################


# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)

# SET(CMAKE_CXX_STANDARD 17)

SET(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_BINARY_DIR}/include)
SET(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

# include_directories(${ROOT_INCLUDE_DIRS}) # I've been typically including this, but it seems unnecessary?
include_directories(
  "${CMAKE_SOURCE_DIR}/display/inc"
  ${DELPHES_INCLUDE_DIR}
)

add_subdirectory(display) # in here, we will call the "ROOT_STANDARD_LIBRARY_PACKAGE" macro, which is provided by ROOT

get_property(__allHeaders GLOBAL PROPERTY ROOT_HEADER_TARGETS)
get_property(__allBuiltins GLOBAL PROPERTY ROOT_BUILTIN_TARGETS)
add_custom_target(move_headers ALL DEPENDS ${__allHeaders} ${__allBuiltins})

