project(GSEA)

####################################################################################################
# Build executable
####################################################################################################

add_executable(gsea main.cpp)
target_link_libraries(gsea gtcore)
set_target_properties(gsea PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(gsea)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS gsea
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
