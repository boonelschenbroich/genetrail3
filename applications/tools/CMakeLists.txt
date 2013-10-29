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

add_executable(cluster_genes cluster_genes.cpp)
target_link_libraries(cluster_genes gtcore gtcluster)
GT2_COMPILE_FLAGS(cluster_genes)

INSTALL(TARGETS mconvert cluster_genes
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
