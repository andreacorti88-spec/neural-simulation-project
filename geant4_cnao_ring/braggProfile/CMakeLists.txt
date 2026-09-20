cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(braggProfile)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(braggProfile braggProfile.cc)
target_link_libraries(braggProfile ${Geant4_LIBRARIES})
