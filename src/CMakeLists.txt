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
option(USING_DOCKER OFF)
if(USING_DOCKER)
	set(SeqAn3_DIR "/deps/seqan3/build_system")
	set(THREAD_POOL_INCLUDE "/deps/ThreadPool")
else()
	set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/deps/seqan3/build_system")
	set(THREAD_POOL_INCLUDE "${CMAKE_SOURCE_DIR}/deps/ThreadPool")
endif(USING_DOCKER)
unset(USING_DOCKER CACHE)

find_package (SeqAn3 REQUIRED)

#################### arg_parse #######################
add_library(arg_parse arg_parse.cpp)
target_include_directories( arg_parse PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(arg_parse LINK_PUBLIC seqan3::seqan3)

#################### file_to_map #######################
add_library(file_to_map file_to_map.cpp)
target_include_directories( file_to_map PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(file_to_map LINK_PUBLIC seqan3::seqan3)


#################### scanner #######################
add_library(scanner fm_scanner.cpp)
target_include_directories( scanner PRIVATE  "${CMAKE_SOURCE_DIR}/include"
                                              "${THREAD_POOL_INCLUDE}")
target_link_libraries(scanner LINK_PUBLIC seqan3::seqan3
                                          arg_parse
                                          file_to_map
                                          )

#################### indexer #######################
add_library(indexer fm_indexer.cpp)
target_include_directories( indexer PRIVATE   "${CMAKE_SOURCE_DIR}/include"
					                                    "${THREAD_POOL_INCLUDE}")
target_link_libraries(indexer LINK_PUBLIC seqan3::seqan3
                                          arg_parse
                                          file_to_map
                                          )

#################### speq  #######################
add_executable (speq main.cpp)
target_include_directories(speq PRIVATE   "${CMAKE_SOURCE_DIR}/include"
					                                "${THREAD_POOL_INCLUDE}")
target_link_libraries(speq LINK_PUBLIC  arg_parse 
                                        indexer
#                                        scanner 
                                        seqan3::seqan3
                                        )

install(TARGETS speq
        DESTINATION bin)
