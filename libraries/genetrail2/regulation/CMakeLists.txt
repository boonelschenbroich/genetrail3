project(GENETRAIL2_REGULATION_LIBRARY)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build GeneTrail2 regulation library.")
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

include_directories(
	"${CMAKE_SOURCE_DIR}/libraries/"
	"${CMAKE_BINARY_DIR}/libraries/"
	"${RapidJSON_INCLUDE_DIR}"
)

add_library(gtregulation SHARED
	${REGULATION_SOURCES}
	${REGULATION_HEADERS}
)

set_target_properties(gtregulation PROPERTIES LINKER_LANGUAGE CXX)

GT2_COMPILE_FLAGS(gtregulation)

target_link_libraries(gtregulation gtcore ${GT2_CORE_DEP_LIBRARIES})

####################################################################################################
# Install the library
####################################################################################################

INSTALL(FILES ${CORE_HEADERS}
	DESTINATION include/genetrail2/regulation/
)

INSTALL(TARGETS gtregulation
	EXPORT GeneTrail2Targets
	RUNTIME DESTINATION bin COMPONENT libraries
	LIBRARY DESTINATION lib COMPONENT libraries
	ARCHIVE DESTINATION lib COMPONENT libraries
)

