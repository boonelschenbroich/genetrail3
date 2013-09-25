project(Tools)

####################################################################################################
# Build executable
####################################################################################################

add_executable(mconvert mconvert.cpp)
target_link_libraries(mconvert gtcore)

INSTALL(TARGETS mconvert
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
