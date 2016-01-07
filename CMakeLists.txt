project(GeneTrail2)

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

####################################################################################################
# Manual entries (edit this for new release)
####################################################################################################

set(GENETRAIL2_PACKAGE_VERSION_MAJOR "0")
set(GENETRAIL2_PACKAGE_VERSION_MINOR "1")
set(GENETRAIL2_PACKAGE_VERSION_PATCH "0")
set(GENETRAIL2_PACKAGE_VERSION "${GENETRAIL2_PACKAGE_VERSION_MAJOR}.${GENETRAIL2_PACKAGE_VERSION_MINOR}.${GENETRAIL2_PACKAGE_VERSION_PATCH}")

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

# Find Boost
find_package(Boost 1.55.0 COMPONENTS iostreams filesystem system program_options REQUIRED)

if(NOT Boost_FOUND)
	message(FATAL_ERROR "Boost version 1.55.0 or higher is required by GeneTrail 2")
endif()

SET(Boost_LIBRARIES
	${Boost_IOSTREAMS_LIBRARY}
	${Boost_FILESYSTEM_LIBRARY}
	${Boost_SYSTEM_LIBRARY}
	${Boost_PROGRAM_OPTIONS_LIBRARY}
)

include_directories(${Boost_INCLUDE_DIRS})

# Find Eigen 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find additional multiprecision libraries
SET(GENETRAIL2_HAS_GMP FALSE)
SET(GENETRAIL2_HAS_MPFR FALSE)

OPTION(USE_GMP "Should the GMP library be used if available" ON)
OPTION(USE_MPFR "Should the MPFR library be used if available" OFF)

if(NOT USE_GMP)
	if(USE_MPFR)
		message(FATAL_ERROR "In order to use MPFR you also need to use GMP")
	endif()
endif()

if(USE_GMP)
	find_package(GMP)
	if(GMP_FOUND)
		include_directories(${GMP_INCLUDE_DIRS})
		set(GENETRAIL2_HAS_GMP TRUE)

		# We use GMP as default. It seems to be faster and we do not
		# need the additional features provided by MPFR.
		if(USE_MPFR)
			find_package(MPFR REQUIRED)
			if(MPFR_FOUND)
				include_directories(${MPFR_INCLUDE_DIRS})
				set(GENETRAIL2_HAS_MPFR TRUE)
			else()
				message(FATAL_ERROR "Could not find the MPFR library")
			endif()
		endif()
	endif()
endif()

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

###################################################################################################
# Add subdirectories
###################################################################################################

add_subdirectory(libraries)
add_subdirectory(applications)

####################################################################################################
# Documentation
####################################################################################################

include(cmake/documentation.cmake)
include(cmake/GT2ConfigHeader.cmake)
include(cmake/GT2Config.cmake)
