####################################################################################################
# Build Tests
####################################################################################################
project(GENETRAIL2_TESTS)


# Configure unit tests. This will only happen if the path to the
# gtest src directory has been set
if(NOT GTEST_SRC_DIR)
    message(STATUS "GTEST_SRC_DIR not set, disabling unit tests")
    return()
endif()

include_directories("${CMAKE_SOURCE_DIR}/libraries")

enable_testing()

add_subdirectory(core)

