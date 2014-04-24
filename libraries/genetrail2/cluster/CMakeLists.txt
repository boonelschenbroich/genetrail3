project(GENETRAIL2_CLUSTER_LIBRARY)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

find_package(METIS)

if(NOT METIS_FOUND)
	message(STATUS "Metis not found, cannot build cluster library")
	return()
endif()

####################################################################################################
# Build a shared library
####################################################################################################

include(files.cmake)

include_directories(
	"${CMAKE_SOURCE_DIR}/libraries/"
	${METIS_INCLUDE_DIRS}
)

add_library(gtcluster SHARED
	${CLUSTER_SOURCES}
	${CLUSTER_HEADERS}
)

GT2_COMPILE_FLAGS(gtcluster)

target_link_libraries(gtcluster
	${METIS_LIBRARIES}
)

####################################################################################################
# Install the library
####################################################################################################

install(FILES ${CLUSTER_HEADERS}
	DESTINATION include/cluster/
)

install(TARGETS gtcluster
	RUNTIME DESTINATION bin COMPONENT libraries
	LIBRARY DESTINATION lib COMPONENT libraries
	ARCHIVE DESTINATION lib COMPONENT libraries
)
