project(Enrich)

####################################################################################################
# Build executable
####################################################################################################

macro(add_enrichment NAME)
	add_executable(${NAME} ${NAME}.cpp)
	target_link_libraries(${NAME} gtcore gtenrichment ${BOOST_LIBRARIES})
	GT2_COMPILE_FLAGS(${NAME})
	set_target_properties(${NAME} PROPERTIES
	    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
	)
endmacro()

add_executable(hotelling_t_test hotelling.cpp)
target_link_libraries(hotelling_t_test gtcore)
set_target_properties(hotelling_t_test PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)
GT2_COMPILE_FLAGS(hotelling_t_test)

add_enrichment(gsea)
add_enrichment(weighted-gsea)

add_enrichment(ora)
configure_file(${CMAKE_SOURCE_DIR}/scripts/enrichment/parallelORA.in ${CMAKE_SOURCE_DIR}/scripts/enrichment/parallelORA @ONLY)
file(COPY ${CMAKE_SOURCE_DIR}/scripts/enrichment/parallelORA
		DESTINATION ${CMAKE_BINARY_DIR}/applications/enrichment/
		FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
		GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

add_enrichment(htests)
add_enrichment(wilcoxon)
add_enrichment(enrichment enrichment.cpp)
add_enrichment(contingency_ora)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS hotelling_t_test gsea ora htests enrichment weighted-gsea contingency_ora
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
