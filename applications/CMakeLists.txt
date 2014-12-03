project(GENETRAIL2_APPLICATIONS)

if(NOT TARGET gtcore)
	return()
endif()

include_directories("${CMAKE_SOURCE_DIR}/libraries/")
include_directories("${CMAKE_BINARY_DIR}/libraries/")

add_subdirectory(enrichment)
add_subdirectory(fidepa)
add_subdirectory(glasso)
add_subdirectory(tools)
add_subdirectory(geo)
add_subdirectory(compute)
