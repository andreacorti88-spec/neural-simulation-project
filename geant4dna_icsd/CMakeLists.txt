# ------------------------------------------------------------------
# nanoICSD -- Geant4-DNA project to reproduce the Ionization Cluster
# Size Distribution (ICSD) test case of:
#   Villagrasa, Baiocco et al., "Evaluation of the uncertainty in
#   calculating nanodosimetric quantities...", PLOS One 21(1):e0340500
#   (2026), DOI 10.1371/journal.pone.0340500
#
# Requires Geant4 (>=10.7, tested design targets 11.x) built WITH
# the G4EMLOW / Geant4-DNA low-energy datasets installed and the
# G4LEDATA (and related G4EM*DATA) environment variables set, as for
# any standard Geant4-DNA example.
# ------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(nanoICSD)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(nanoICSD nanoICSD.cc ${sources} ${headers})
target_link_libraries(nanoICSD ${Geant4_LIBRARIES})

file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS nanoICSD DESTINATION bin)
