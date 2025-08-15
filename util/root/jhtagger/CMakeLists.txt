# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.5...3.30)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
  jhtagger
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

# Find FastJet.
# Will search CMAKE_PREFIX_PATH, plus the FASTJET_EXT_DIR
# that we define here (to pick up a possible local installation).
set(FASTJET_EXT_DIR "${REPO_DIR}/external/fastjet/fastjet-install")

set(FASTJET_INCLUDE_DIR "${FASTJET_EXT_DIR}/include")
if(EXISTS "${FASTJET_INCLUDE_DIR}/fastjet/ClusterSequence.hh")
    message(STATUS "Found FastJet include directory: ${FASTJET_INCLUDE_DIR}")
else()
    set(FASTJET_INCLUDE_DIR "")
    message(WARNING "FastJet headers not found in expected location")
endif()

find_library(FASTJET_LIB
    NAMES fastjet libfastjet
    PATHS
        ${FASTJET_EXT_DIR}/lib
        ${FASTJET_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT FASTJET_LIB)
  find_library(FASTJET_LIB
      NAMES fastjet libfastjet
      PATHS
          ${FASTJET_EXT_DIR}/lib
          ${FASTJET_EXT_DIR}
  )
endif()

message(STATUS "Found FastJet library: ${FASTJET_LIB}")

# We also need the fastjet tools library.
find_library(FASTJET_TOOLS_LIB
    NAMES fastjettools libfastjettools
    PATHS
        ${FASTJET_EXT_DIR}/lib
        ${FASTJET_EXT_DIR}
    NO_DEFAULT_PATH
)

if(NOT FASTJET_TOOLS_LIB)
  find_library(FASTJET_TOOLS_LIB
      NAMES fastjettools libfastjettools
      PATHS
          ${FASTJET_EXT_DIR}/lib
          ${FASTJET_EXT_DIR}
  )
endif()

message(STATUS "Found FastJet tools library: ${FASTJET_TOOLS_LIB}")

# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)


SET(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_BINARY_DIR}/include)
SET(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

# include_directories(${ROOT_INCLUDE_DIRS})
include_directories(
  "${CMAKE_SOURCE_DIR}/jhtagger/inc"
  ${FASTJET_INCLUDE_DIR}
)

# # Add the Fastjet include libraries.
# include_directories($ENV{FASTJET_INCLUDE_DIR}) # This environment variable has been created externally, by the (shell) build script.
# message(STATUS "Adding Fastjet include directory: $ENV{FASTJET_INCLUDE_DIR}") # TODO: Could make this fancier, give a different message if the find actually fails.

add_subdirectory(jhtagger)

get_property(__allHeaders GLOBAL PROPERTY ROOT_HEADER_TARGETS)
get_property(__allBuiltins GLOBAL PROPERTY ROOT_BUILTIN_TARGETS)
add_custom_target(move_headers ALL DEPENDS ${__allHeaders} ${__allBuiltins})
