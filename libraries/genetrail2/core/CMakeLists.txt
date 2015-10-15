project(GENETRAIL2_CORE_LIBRARY)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build GeneTrail2 core library.")
	return()
endif()

if(NOT RapidJSON_INCLUDE_DIR)
	message(STATUS "RapidJSON_INCLUDE_DIR not specfied, cannot build GeneTrail2 core library.")
	return()
endif()

SET(GT2_CORE_DEP_LIBRARIES
	${Boost_LIBRARIES}
)

if(GENETRAIL2_HAS_GMP)
	list(APPEND GT2_CORE_DEP_LIBRARIES ${GMP_LIBRARIES})
endif()

if(GENETRAIL2_HAS_MPFR)
	list(APPEND GT2_CORE_DEP_LIBRARIES ${MPFR_LIBRARIES})
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

target_include_directories(gtcore PRIVATE "${CMAKE_BINARY_DIR}/libraries/" "${RapidJSON_INCLUDE_DIR}")

target_link_libraries(gtcore ${GT2_CORE_DEP_LIBRARIES})

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

