# This CMakeLists file runs a ROOT CMake macro -- this has been set up by running "find_package (ROOT ...)" in the main CMakeLists file.
# (This functionality is available at least since ROOT v. 6.22, if not earlier).

ROOT_STANDARD_LIBRARY_PACKAGE(EventDisplay
  HEADERS
    display/DisplayInterface.h
    display/EventDisplay.h
    display/ExDelphesDisplay.h
    display/BranchElement.h
    display/CaloData.h
    display/TestUtils.h
  SOURCES
    src/DisplayInterface.cxx
    src/EventDisplay.cxx
    src/ExDelphesDisplay.cxx
    src/BranchElement.cxx
    src/CaloData.cxx
    src/TestUtils.cxx
  LINKDEF
    Display_LinkDef.h
  DEPENDENCIES
    ROOT::Core ROOT::Tree ROOT::EG ROOT::Geom ROOT::Eve ROOT::Gui ROOT::GuiHtml
  DICTIONARY_OPTIONS
    -writeEmptyRootPCM
)

# TODO: Not sure if this is really the best way to go,
# the above ROOT CMake macro required include_directories()
# to be called on ${CMAKE_SOURCE_DIR}/display/inc

# Add the include directories.
message(STATUS "Adding to target_include_directories: ${DELPHES_INCLUDE_DIR}")
target_include_directories(EventDisplay
    PRIVATE
        "${CMAKE_SOURCE_DIR}/display/inc"
        ${DELPHES_INCLUDE_DIR}
)
# Add the Delphes libraries that the ROOT macro doesn't handle
target_link_libraries(EventDisplay
    PRIVATE
        ${DELPHES_LIBRARIES}
        ${DELPHES_DISPLAY_LIBRARIES}
        ROOT::Core
        ROOT::Tree
        ROOT::EG
        ROOT::Geom
        ROOT::Eve
        ROOT::Gui
        ROOT::GuiHtml
        # ROOT::GeomPainter
        # ROOT::GeomBuilder
)
