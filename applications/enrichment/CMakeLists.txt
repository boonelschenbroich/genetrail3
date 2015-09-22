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

add_executable(hotelling_t_test main.cpp)
target_link_libraries(hotelling_t_test gtcore)
set_target_properties(hotelling_t_test PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)
GT2_COMPILE_FLAGS(hotelling_t_test)

add_enrichment(gsea)
add_enrichment(weighted-gsea)
add_enrichment(ora)
add_enrichment(htests)
add_enrichment(enrichment enrichment.cpp)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS hotelling_t_test gsea ora htests enrichment weighted-gsea
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
