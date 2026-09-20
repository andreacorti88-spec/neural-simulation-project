# ------------------------------------------------------------------
# cnaoRingImpact -- physical dose/hit map of a CNAO-style clinical
# ion beam (proton or carbon-12) on a ring of 100 tissue-equivalent
# "neuron" spheres, the physical counterpart of the N_E=100
# excitatory ring in the spiking ring attractor of neuroni-progetto.
#
# Uses QBBC (standard Geant4 physics for clinical proton/ion
# energies), NOT the Geant4-DNA track-structure physics of the
# sibling project nanoICSD (valid only up to ~1 MeV). Requires a
# standard Geant4 11.x build; no Geant4-DNA/G4EMLOW dataset needed.
# ------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(cnaoRingImpact)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(cnaoRingImpact cnaoRingImpact.cc ${sources} ${headers})
target_link_libraries(cnaoRingImpact ${Geant4_LIBRARIES})

file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS cnaoRingImpact DESTINATION bin)
