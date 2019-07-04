cmake_minimum_required (VERSION 2.8)
cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

add_definitions(--std=c++11 )

project (WaveSolver CXX)

set( LEARN_CMAKE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set( LEARN_CMAKE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} )

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/data/")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data")
    message("Created directory ${CMAKE_CURRENT_BINARY_DIR}/data")
else()
    # file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/data/") #doesn't work somehow
    file(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/data/*")
    foreach(file ${files})
        file(REMOVE "${file}")
    endforeach()
    message("Cleaned previous output from directory ${CMAKE_CURRENT_BINARY_DIR}/data/*")
endif()

include_directories ("${PROJECT_SOURCE_DIR}/../")
include_directories ("${PROJECT_SOURCE_DIR}/Catch2")
include_directories ("${PROJECT_SOURCE_DIR}/Acoustic2d")
include_directories ("${PROJECT_SOURCE_DIR}/Test")


add_library(Acoustic2d  Acoustic2d/Acoustic2d.cpp)

#set (LEARN_CMAKE_VERSION_MAJOR 1)
#set (LEARN_CMAKE_VERSION_MINOR 0)
#set (LEARN_CMAKE_VERSION_PATCH 0)

add_executable(solver main.cpp)
target_link_libraries (solver Acoustic2d)
