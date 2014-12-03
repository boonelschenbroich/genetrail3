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

include(files.cmake)

add_library(gtcore SHARED
	${CORE_SOURCES}
	${CORE_HEADERS}
)

GT2_COMPILE_FLAGS(gtcore)

target_include_directories(gtcore PRIVATE "${CMAKE_BINARY_DIR}/libraries/")
target_link_libraries(gtcore
	${Boost_LIBRARIES}
)

####################################################################################################
# Install the library
####################################################################################################

INSTALL(FILES ${CORE_HEADERS}
	DESTINATION include/genetrail2/core/
)

INSTALL(TARGETS gtcore
	EXPORT GeneTrail2Targets
	RUNTIME DESTINATION bin COMPONENT libraries
	LIBRARY DESTINATION lib COMPONENT libraries
	ARCHIVE DESTINATION lib COMPONENT libraries
)

