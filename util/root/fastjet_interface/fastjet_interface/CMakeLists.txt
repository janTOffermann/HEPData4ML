# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(FastjetInterface
  HEADERS
    fastjet_interface/FastjetInterface.h
  SOURCES
    src/FastjetInterface.cxx
  LINKDEF
    FastjetInterface_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

target_include_directories(FastjetInterface PRIVATE
  "${FASTJET_INCLUDE_DIR}"
    "${CMAKE_SOURCE_DIR}/fastjet_interface/inc"

)
target_link_libraries(FastjetInterface PRIVATE
  "${FASTJET_LIB}"
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
)
