project(GENETRAIL2_ENRICHMENT_LIBRARY)

SET(GT2_ENRICHMNET_DEP_LIBRARIES
	${Boost_LIBRARIES}
)

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "ENRICHMENT")

add_header_to_library(
	PermutationTest.h
)

add_to_library(
	common
	CommandLineInterface
	EnrichmentAlgorithm
	Parameters
	SetLevelStatistics
)

add_library(gtenrichment SHARED
	${ENRICHMENT_SOURCES}
	${ENRICHMENT_HEADERS}
)

GT2_COMPILE_FLAGS(gtenrichment)

target_include_directories(gtenrichment PRIVATE
	"${CMAKE_BINARY_DIR}/libraries/"
	"${CMAKE_SOURCE_DIR}/libraries/"
)

target_link_libraries(gtenrichment gtcore ${GT2_ENRICHMENT_DEP_LIBRARIES})

####################################################################################################
# Install the library
####################################################################################################

INSTALL(FILES ${ENRICHMENT_HEADERS}
	DESTINATION include/genetrail2/enrichment/
)

INSTALL(TARGETS gtenrichment
	EXPORT GeneTrail2Targets
	RUNTIME DESTINATION bin COMPONENT libraries
	LIBRARY DESTINATION lib COMPONENT libraries
	ARCHIVE DESTINATION lib COMPONENT libraries
)

