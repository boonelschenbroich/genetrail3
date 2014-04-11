project(GSEA)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build wilcoxon")
    return()
endif()

####################################################################################################
# Build executable
####################################################################################################

add_executable(wilcoxon main.cpp)
target_link_libraries(wilcoxon gtcore)
set_target_properties(wilcoxon PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(wilcoxon)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS wilcoxon
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
