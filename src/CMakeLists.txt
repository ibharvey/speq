cmake_minimum_required (VERSION 3.5)

project (speq CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

################## Add external projects #####################
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

################## Add zlib #####################
find_package (ZLIB REQUIRED)

################## Add seqan3 #####################
set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/deps/seqan3/build_system")
find_package (SeqAn3 REQUIRED)


#################### arg_parse #######################
add_library(arg_parse arg_parse.cpp)
target_include_directories( arg_parse PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(arg_parse LINK_PUBLIC seqan3::seqan3)

#################### file_to_map #######################
add_library(file_to_map file_to_map.cpp)
target_include_directories( file_to_map PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(file_to_map LINK_PUBLIC seqan3::seqan3)


#################### aligner #######################
add_library(aligner aligner.cpp)
target_include_directories( aligner PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(aligner LINK_PUBLIC seqan3::seqan3
                                          arg_parse
                                          file_to_map
                                          )

#################### speq  #######################
add_executable (speq main.cpp)
target_include_directories(speq PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(speq LINK_PUBLIC  arg_parse 
                                        aligner 
                                        seqan3::seqan3
                                        )

install(TARGETS speq
        DESTINATION bin)