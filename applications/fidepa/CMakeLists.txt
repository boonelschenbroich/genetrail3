project(FiDePa)

####################################################################################################
# Build executable
####################################################################################################

add_executable(fidepa fidepa.cpp)
target_link_libraries(fidepa gtcore)
set_target_properties(fidepa PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(fidepa)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS fidepa
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
