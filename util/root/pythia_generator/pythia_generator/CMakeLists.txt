# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(PythiaGenerator
  HEADERS
    pythia_generator/PythiaGenerator.h
  SOURCES
    src/PythiaGenerator.cxx
  LINKDEF
    PythiaGenerator_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

target_include_directories(PythiaGenerator PRIVATE
  "${PYTHIA8_INCLUDE_DIR}"
    "${CMAKE_SOURCE_DIR}/pythia_generator/inc"

)
target_link_libraries(PythiaGenerator PRIVATE
  "${PYTHIA8_LIB}"
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
)
