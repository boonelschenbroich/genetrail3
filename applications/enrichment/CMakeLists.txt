project(Enrich)

####################################################################################################
# Build executable
####################################################################################################

add_executable(enrich main.cpp)
target_link_libraries(enrich gtcore)
set_target_properties(enrich PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)
GT2_COMPILE_FLAGS(enrich)

add_library(enrichment_common STATIC common.cpp common.h)
GT2_COMPILE_FLAGS(enrichment_common)

add_executable(gsea gsea.cpp)
target_link_libraries(gsea gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(gsea)
set_target_properties(gsea PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(ora ora.cpp)
target_link_libraries(ora gtcore enrichment_common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(ora)
set_target_properties(ora PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS enrich gsea ora
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
