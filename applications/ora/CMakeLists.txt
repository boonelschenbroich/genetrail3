project(ORA)

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

add_executable(ora main.cpp)
target_link_libraries(ora gtcore)
set_target_properties(ora PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(ora)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS ora
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
