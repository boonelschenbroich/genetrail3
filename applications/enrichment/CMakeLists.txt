project(Enrich)

####################################################################################################
# Build executable
####################################################################################################

add_executable(hotelling_t_test main.cpp)
target_link_libraries(hotelling_t_test gtcore)
set_target_properties(hotelling_t_test PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)
GT2_COMPILE_FLAGS(hotelling_t_test)

add_library(enrichment_common STATIC common.cpp common.h)
GT2_COMPILE_FLAGS(enrichment_common)
set_target_properties(enrichment_common PROPERTIES
	INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(gsea gsea.cpp)
target_link_libraries(gsea gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(gsea)
set_target_properties(gsea PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(weighted-gsea weighted-gsea.cpp)
target_link_libraries(weighted-gsea gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(weighted-gsea)
set_target_properties(weighted-gsea PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(ora ora.cpp)
target_link_libraries(ora gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(ora)
set_target_properties(ora PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(htests htests.cpp)
target_link_libraries(htests gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(htests)
set_target_properties(htests PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(enrichment enrichment.cpp)
target_link_libraries(enrichment gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(enrichment)
set_target_properties(enrichment PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS hotelling_t_test gsea ora htests enrichment weighted-gsea
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
