project(Enrich)

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

add_executable(enrich main.cpp)
target_link_libraries(enrich gtcore)
set_target_properties(enrich PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(enrich)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS enrich
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)