project(GENETRAIL2_CORE_LIBRARY)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build GeneTrail2 core library")
    return()
endif()

####################################################################################################
# Build a shared library
####################################################################################################

include(src/files.cmake)

add_library(gtcore SHARED 
	${CORE_SOURCES} 
	${CORE_HEADERS}
)

target_link_libraries(gtcore
	${Boost_LIBRARIES}
)

####################################################################################################
# Install the library
####################################################################################################

INSTALL(FILES ${CORE_HEADERS}
    DESTINATION include/core/
)

INSTALL(TARGETS gtcore
    RUNTIME DESTINATION bin COMPONENT libraries
    LIBRARY DESTINATION lib COMPONENT libraries
    ARCHIVE DESTINATION lib COMPONENT libraries
)

####################################################################################################
# Build Tests
####################################################################################################

# Configure unit tests. This will only happen if the path to the
# gtest src directory has been set

if(NOT GTEST_SRC_DIR)
    message(STATUS "GTEST_SRC_DIR not set, disabling unit tests")
	return()
endif()

enable_testing()

add_subdirectory(test)

