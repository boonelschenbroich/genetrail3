project(Tools)

####################################################################################################
# Build executable
####################################################################################################

include_directories(
	"${CMAKE_SOURCE_DIR}/libraries/core/src/"
	"${CMAKE_SOURCE_DIR}/libraries/cluster/src/"
)

add_executable(mconvert mconvert.cpp)
target_link_libraries(mconvert gtcore)
GT2_COMPILE_FLAGS(mconvert)
list(APPEND EXECUTABLES mconvert)

add_executable(filter-scores filter_scores.cpp)
target_link_libraries(filter-scores gtcore)
GT2_COMPILE_FLAGS(filter-scores)
list(APPEND EXECUTABLES filter-scores)

if(TARGET gtcluster)
	add_executable(cluster_genes cluster_genes.cpp)
	target_link_libraries(cluster_genes gtcore gtcluster)
	GT2_COMPILE_FLAGS(cluster_genes)
	list(APPEND EXECUTABLES cluster_genes)
endif()

INSTALL(TARGETS ${EXECUTABLES}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
