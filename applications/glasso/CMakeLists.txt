project(Glasso)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
    message(STATUS "Boost not found, cannot build FiDePa")
    return()
endif()

include_directories(
	"${CMAKE_SOURCE_DIR}/libraries/core/src/"
)

####################################################################################################
# Build executable
####################################################################################################

enable_language(Fortran)

add_library(lglasso STATIC glasso.f90)
set_target_properties(lglasso PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(glasso main.cpp glasso.h)
target_link_libraries(glasso gtcore lglasso)
set_target_properties(glasso PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(glasso)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS glasso
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
