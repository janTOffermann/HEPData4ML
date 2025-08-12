# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(JHTagger
  HEADERS
    jhtagger/JHTagger.h
  SOURCES
    src/JHTagger.cxx
  LINKDEF
    JHTagger_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::RIO ROOT::MathCore
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

target_include_directories(JHTagger PRIVATE
  "${FASTJET_INCLUDE_DIR}"
    "${CMAKE_SOURCE_DIR}/jhtagger/inc"

)
target_link_libraries(JHTagger PRIVATE
  "${FASTJET_LIB}"
  "${FASTJET_TOOLS_LIB}"
  ROOT::Core
  ROOT::RIO
  ROOT::MathCore
)
