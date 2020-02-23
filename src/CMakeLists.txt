cmake_minimum_required (VERSION 3.4)

project (hopa CXX)


# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/deps/seqan3/build_system")
find_package (SeqAn3 REQUIRED)

add_executable (hopa main.cpp)


# arg_parse library
add_library(arg_parse arg_parse.cpp)
target_include_directories( arg_parse PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )
target_link_libraries(arg_parse LINK_PUBLIC seqan3::seqan3)

# traversal library
add_library(traversal traversal.cpp)

target_include_directories( traversal PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )
target_link_libraries(traversal LINK_PUBLIC seqan3::seqan3)

target_include_directories(hopa PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
  )

target_link_libraries(hopa LINK_PUBLIC arg_parse traversal seqan3::seqan3)

