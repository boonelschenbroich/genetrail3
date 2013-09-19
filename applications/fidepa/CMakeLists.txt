project(FiDePa)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
    MESSAGE(STATUS "Boost not found, cannot build FiDePa")
    return()
endif()

####################################################################################################
# Build executable
####################################################################################################

add_executable(fidepa fidepa.cpp)
target_link_libraries(fidepa gtcore)
set_target_properties(fidepa PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Build executable
####################################################################################################

INSTALL(TARGETS fidepa
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
