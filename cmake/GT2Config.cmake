list(APPEND GENETRAIL2_LIBRARY_INCLUDE_DIRS
	${EIGEN3_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
	"${PROJECT_SOURCE_DIR}/libraries"
	"${PROJECT_BINARY_DIR}/libraries"
)

# TODO: This can probably be changed when we switch to CMake 3.0.
#       That version offers to automatically detect the built targets
#       in the export directive.
if(TARGET gtcore)
	list(APPEND GENETRAIL2_AVAILABLE_LIBRARIES gtcore)
endif()

if(TARGET gtcluster)
	list(APPEND GENETRAIL2_AVAILABLE_LIBRARIES gtcluster)
endif()

export(TARGETS ${GENETRAIL2_AVAILABLE_LIBRARIES}
	FILE "${PROJECT_BINARY_DIR}/GeneTrail2Targets.cmake"
)

export(PACKAGE GeneTrail2)

configure_file(cmake/GeneTrail2Config.cmake.in "${PROJECT_BINARY_DIR}/GeneTrail2Config.cmake" @ONLY)
configure_file(cmake/GeneTrail2ConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/GeneTrail2ConfigVersion.cmake" @ONLY)
