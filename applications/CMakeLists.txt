project(GENETRAIL2_APPLICATIONS)

include_directories("${CMAKE_SOURCE_DIR}/libraries/")

add_subdirectory(enrichment)
add_subdirectory(fidepa)
add_subdirectory(glasso)
add_subdirectory(ora)
add_subdirectory(tools)
add_subdirectory(gsea)
