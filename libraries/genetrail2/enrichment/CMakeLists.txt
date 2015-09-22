project(GENETRAIL2_ENRICHMENT_LIBRARY)

SET(GT2_ENRICHMNET_DEP_LIBRARIES
	${Boost_LIBRARIES}
)

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "ENRICHMENT")

add_to_library(
	common
	CommandLineInterface
	Parameters
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
