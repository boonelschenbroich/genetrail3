project(GENETRAIL2)

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

####################################################################################################
# Manual entries (edit this for new release)
####################################################################################################

set(GENETRAIL2_PACKAGE_VERSION_MAJOR "2")
set(GENETRAIL2_PACKAGE_VERSION_MINOR "0")
set(GENETRAIL2_PACKAGE_VERSION_PATCH "0")

####################################################################################################
# Set the path to the additional cmake files.
####################################################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

include(GT2CompilerSpecific)
include(GT2Macros)

####################################################################################################
# Build Tests
####################################################################################################

# This call needs to be in the top level file before any
# add_subdirectory calls in order to allow the testing
# framework to pick up all unit test.

# We require the gtest src directory to be set
if(NOT GTEST_SRC_DIR)
    message(STATUS "GTEST_SRC_DIR not set, disabling unit tests")
	SET(GT2_ENABLE_TESTING FALSE)
else()
	enable_testing()
	SET(GT2_ENABLE_TESTING TRUE)
endif()

####################################################################################################
# Find third party libraries
####################################################################################################

find_package(Boost 1.55.0 COMPONENTS iostreams regex filesystem system program_options REQUIRED)

if(NOT Boost_FOUND)
	message(FATAL_ERROR "Boost version 1.55.0 or higher is required by GeneTrail 2")
endif()

SET(Boost_LIBRARIES
	${Boost_IOSTREAMS_LIBRARY}
	${Boost_REGEX_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

####################################################################################################
# Check if a build type is given
####################################################################################################

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

####################################################################################################
# Check if the rpath should be used to install GeneTrail2
####################################################################################################

option(USE_INSTALL_RPATH "Should CMake set the install rpath?" OFF)

if(USE_INSTALL_RPATH)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

####################################################################################################
# Compile Flags
####################################################################################################

# TODO set visibility = hidden as standard 
# There is a Makro in the BALL project to do this

###################################################################################################
# Build libraries
###################################################################################################

add_subdirectory(libraries)

###################################################################################################
# Build Applications
###################################################################################################

add_subdirectory(applications)

add_subdirectory(bindings)

####################################################################################################
# Documentation
####################################################################################################

include(cmake/documentation.cmake)


