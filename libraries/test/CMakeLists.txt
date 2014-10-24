####################################################################################################
# Build Tests
####################################################################################################
project(GENETRAIL2_TESTS)

add_library(test_driver main_test.cpp)
include(cmake/GeneTrail2TestExecuter.cmake)

include_directories("${GTEST_SRC_DIR}/include" "${CMAKE_SOURCE_DIR}/libraries")

add_subdirectory(core)
