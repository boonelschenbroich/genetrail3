project(GSEA)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build Enrich")
    return()
endif()

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
