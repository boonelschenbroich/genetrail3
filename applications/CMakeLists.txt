project(GENETRAIL2_APPLICATIONS)

if(NOT TARGET gtcore)
	return()
endif()

add_library(applicationCommon STATIC matrixTools.h matrixTools.cpp)
target_include_directories(applicationCommon PRIVATE
	${CMAKE_SOURCE_DIR}/libraries/genetrail2/core
	${CMAKE_BINARY_DIR}/libraries/genetrail2/core
)
GT2_COMPILE_FLAGS(applicationCommon)

include_directories("${CMAKE_SOURCE_DIR}/libraries/")
include_directories("${CMAKE_BINARY_DIR}/libraries/")

add_subdirectory(enrichment)
add_subdirectory(fidepa)
add_subdirectory(glasso)
add_subdirectory(tools)
add_subdirectory(geo)
add_subdirectory(compute)
add_subdirectory(test)
