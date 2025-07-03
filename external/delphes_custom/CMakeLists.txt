cmake_minimum_required(VERSION 3.10)

project(DelphesHepMC3ROOT)

# Set C++ standard to match ROOT/Delphes
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compile options similar to original Delphes
add_compile_options(-fPIC -Wall -D_FILE_OFFSET_BITS=64)

# Define the path to this repository -- relative to this CMake file.
# NOTE: If this executable is made into a separate git submodule,
#       some of this will likely have to change.
#       Maybe have external installations modify CMAKE_PREFIX_PATH?
set(REPO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")

#######################################################
#########
# ROOT
#########

# Find ROOT - required for Delphes
find_package(ROOT CONFIG REQUIRED)
#######################################################
#########
# HEPMC3
#########

# Find HepMC3 -- required for the HepMC3/ROOT format handling.
# Will search CMAKE_PREFIX_PATH, plus the HEPMC3_EXT_DIR
# that we define here (to pick up a possible local installation).
# Note that we're also searching for the HepMC3/ROOT interface.
set(HEPMC3_EXT_DIR "${REPO_DIR}/external/hepmc3")

find_path(HEPMC3_INCLUDE_DIR
    NAMES HepMC3
    PATHS
        ${HEPMC3_EXT_DIR} # fallback
)

find_library(HEPMC3_LIBRARIES
    NAMES HepMC3 libHepMC3
    PATHS
        ${HEPMC3_EXT_DIR}/lib64
        ${HEPMC3_EXT_DIR}/lib
        ${HEPMC3_EXT_DIR}
        ${CMAKE_PREFIX_PATH}
        ${CMAKE_PREFIX_PATH}/lib
        ${CMAKE_PREFIX_PATH}/lib64
)

find_library(HEPMC3_ROOT_LIBRARIES
    NAMES HepMC3rootIO libHepMC3rootIO
    PATHS
        ${HEPMC3_EXT_DIR}/lib64
        ${HEPMC3_EXT_DIR}
        ${CMAKE_PREFIX_PATH}
        ${CMAKE_PREFIX_PATH}/lib
        ${CMAKE_PREFIX_PATH}/lib64
)

message(STATUS "Found HepMC3 headers: ${HEPMC3_INCLUDE_DIR}")
message(STATUS "Found HepMC3 library: ${HEPMC3_LIBRARIES}")
message(STATUS "Found HepMC3/ROOT library: ${HEPMC3_ROOT_LIBRARIES}")

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
        message(FATAL_ERROR "DelphesClasses.h not found in ${DELPHES_EXT_DIR} or ${DELPHES_EXT_DIR}/include")
    endif()
endif()

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

# Check if Delphes was found
if(NOT DELPHES_INCLUDE_DIR)
    message(FATAL_ERROR "Delphes headers not found. Please set DELPHES_EXT_DIR environment variable or CMAKE_PREFIX_PATH")
endif()

if(NOT DELPHES_LIBRARIES)
    message(FATAL_ERROR "Delphes library not found. Please set DELPHES_EXT_DIR environment variable or CMAKE_PREFIX_PATH")
endif()

message(STATUS "Found Delphes headers: ${DELPHES_INCLUDE_DIR}")
message(STATUS "Found Delphes library: ${DELPHES_LIBRARIES}")

#######################################################

# Create the executable with all source files
add_executable(DelphesHepMC3ROOT
    src/DelphesHepMC3ROOT.cpp
    src/DelphesHepMC3ROOTReader.cpp
)

target_include_directories(DelphesHepMC3ROOT
    PRIVATE
        "${CMAKE_SOURCE_DIR}/inc"
        ${DELPHES_INCLUDE_DIR}
        ${HEPMC3_INCLUDE_DIR}
)

target_link_libraries(DelphesHepMC3ROOT
    PRIVATE
        ${DELPHES_LIBRARIES}
        ${HEPMC3_LIBRARIES}
        ${HEPMC3_ROOT_LIBRARIES}
        ROOT::Core
        ROOT::Tree
        ROOT::EG
)

# # Installation
# install(TARGETS DelphesHepMC3ROOT DESTINATION bin)
