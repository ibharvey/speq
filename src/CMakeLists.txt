cmake_minimum_required (VERSION 3.5)

project (speq CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/deps/seqan3/build_system")
find_package (SeqAn3 REQUIRED)
find_package (ZLIB REQUIRED)

add_executable (speq main.cpp)


# arg_parse library
add_library(arg_parse arg_parse.cpp)
target_include_directories( arg_parse PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )
target_link_libraries(arg_parse LINK_PUBLIC seqan3::seqan3)

# aligner library
add_library(aligner aligner.cpp)

target_include_directories( aligner PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )
target_link_libraries(aligner LINK_PUBLIC seqan3::seqan3)

target_include_directories(speq PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )

target_link_libraries(speq LINK_PUBLIC arg_parse aligner seqan3::seqan3)

