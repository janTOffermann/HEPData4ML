# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.1...3.3)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(
  vectorcalcs
  VERSION 1.0
  LANGUAGES C CXX)

# Fetch ROOT! As of ~6.20, this also sets up our ROOT CMake macros.
find_package(ROOT ROOT_VERSION CONFIG REQUIRED)

# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)

SET(CMAKE_CXX_STANDARD 17)

SET(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_BINARY_DIR}/include)
SET(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

include_directories(${ROOT_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/vectorcalcs/inc")

add_subdirectory(vectorcalcs)

get_property(__allHeaders GLOBAL PROPERTY ROOT_HEADER_TARGETS)
get_property(__allBuiltins GLOBAL PROPERTY ROOT_BUILTIN_TARGETS)
add_custom_target(move_headers ALL DEPENDS ${__allHeaders} ${__allBuiltins})
