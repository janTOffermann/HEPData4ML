# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.5...3.30)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
  pileup
  VERSION 1.0
  LANGUAGES C CXX)

# Define the path to this repository -- relative to this CMake file.
# NOTE: If this executable is made into a separate git submodule,
#       some of this will likely have to change.
#       Maybe have external installations modify CMAKE_PREFIX_PATH?
set(REPO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../..")

# Fetch ROOT! As of ~6.20, this also sets up our ROOT CMake macros.
find_package(ROOT CONFIG REQUIRED)

# Find the root-config utility
find_program(ROOT_CONFIG_CMD root-config)

if(ROOT_CONFIG_CMD)
    # Get ROOT's compile flags which include the C++ standard
    execute_process(
        COMMAND ${ROOT_CONFIG_CMD} --cflags
        OUTPUT_VARIABLE ROOT_CFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE ROOT_CONFIG_RESULT
    )

    if(ROOT_CONFIG_RESULT EQUAL 0)
        # message(STATUS "ROOT CFLAGS: ${ROOT_CFLAGS}")

        # Parse for -std=c++XX in the CXXFLAGS
        if(ROOT_CFLAGS MATCHES "-std=c\\+\\+([0-9]+)")
            set(ROOT_DETECTED_STD ${CMAKE_MATCH_1})
            # message(STATUS "Detected C++ standard from ROOT: C++${ROOT_DETECTED_STD}")
            set(CMAKE_CXX_STANDARD ${ROOT_DETECTED_STD})
        else()
            # Fallback if no -std flag found
            set(CMAKE_CXX_STANDARD 17)
            message(WARNING "Could not parse C++ standard from ROOT flags, using C++17")
        endif()
    else()
        message(WARNING "Failed to run root-config --cflags")
        set(CMAKE_CXX_STANDARD 17)
    endif()
else()
    message(WARNING "root-config not found")
    set(CMAKE_CXX_STANDARD 17)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)

SET(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_BINARY_DIR}/include)
SET(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

#########
# HEPMC3
#########

# Find HepMC3 -- required for the HepMC3/ROOT format handling.
# Will first search HEPMC3_EXT_DIR, plus the CMAKE_PREFIX_PATH.
# NOTE: It's important that we first search HEPMC3_EXT_DIR, because
# we'll be loading this package with (Py)ROOT and we might run
# into issues if its built against a different HepMC3 than the
# one we're using in Python.
# Note that we're also searching for the HepMC3/ROOT interface.
set(HEPMC3_EXT_DIR "${REPO_DIR}/external/hepmc/hepmc3-install")
message(STATUS "Set HEPMC3_EXT_DIR = ${HEPMC3_EXT_DIR}")

find_path(HEPMC3_INCLUDE_DIR
    NAMES HepMC3
    PATHS
        ${HEPMC3_EXT_DIR}
        ${HEPMC3_EXT_DIR}/include
    NO_DEFAULT_PATH
)
if(NOT HEPMC3_INCLUDE_DIR)
    find_path(HEPMC3_INCLUDE_DIR
        NAMES HepMC3
    )
endif()


find_library(HEPMC3_LIBRARIES
    NAMES HepMC3 libHepMC3
    PATHS
        ${HEPMC3_EXT_DIR}/lib64
        ${HEPMC3_EXT_DIR}/lib
        ${HEPMC3_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT HEPMC3_LIBRARIES)
    find_library(HEPMC3_LIBRARIES
        NAMES HepMC3 libHepMC3
        PATHS
            ${CMAKE_PREFIX_PATH}
            ${CMAKE_PREFIX_PATH}/lib
            ${CMAKE_PREFIX_PATH}/lib64
    )
endif()

find_library(HEPMC3_ROOT_LIBRARIES
    NAMES HepMC3rootIO libHepMC3rootIO
    PATHS
        ${HEPMC3_EXT_DIR}/lib64
        ${HEPMC3_EXT_DIR}/lib
        ${HEPMC3_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT HEPMC3_ROOT_LIBRARIES)
    find_library(HEPMC3_ROOT_LIBRARIES
        NAMES HepMC3rootIO libHepMC3rootIO
        PATHS
            ${CMAKE_PREFIX_PATH}
            ${CMAKE_PREFIX_PATH}/lib
            ${CMAKE_PREFIX_PATH}/lib64
    )
endif()

message(STATUS "Found HepMC3 headers: ${HEPMC3_INCLUDE_DIR}")
message(STATUS "Found HepMC3 library: ${HEPMC3_LIBRARIES}")
message(STATUS "Found HepMC3/ROOT library: ${HEPMC3_ROOT_LIBRARIES}")

include_directories(
  "${CMAKE_SOURCE_DIR}/pileup/inc"
  ${HEPMC3_INCLUDE_DIR}
)

add_subdirectory(pileup)

get_property(__allHeaders GLOBAL PROPERTY ROOT_HEADER_TARGETS)
get_property(__allBuiltins GLOBAL PROPERTY ROOT_BUILTIN_TARGETS)
add_custom_target(move_headers ALL DEPENDS ${__allHeaders} ${__allBuiltins})
