cmake_minimum_required(VERSION 3.5)
project(speq VERSION 0.1.0)

include(CTest)
enable_testing()

############### Add source ##################
add_subdirectory(src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
