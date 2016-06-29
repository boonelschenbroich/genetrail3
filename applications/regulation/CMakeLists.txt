project(Regulation)

####################################################################################################
# Build executable
####################################################################################################

add_executable(reggae reggae.cpp)
target_link_libraries(reggae applicationCommon gtcore gtregulation)
set_target_properties(reggae PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(reggae)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS reggae
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
