project(GENETRAIL2)

cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)

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

####################################################################################################
# Find third party libraries
####################################################################################################

find_package(Boost 1.46.1 COMPONENTS iostreams regex filesystem system program_options REQUIRED)
find_package(CLN)
find_package(CPPUnit)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CLN_INCLUDE_DIR})

link_directories(${Boost_LIBRARY_DIRS})

include(GT2CompilerSpecific)

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

####################################################################################################
# Documentation
####################################################################################################

include(cmake/documentation.cmake)
