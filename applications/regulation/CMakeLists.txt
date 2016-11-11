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

add_executable(reggae_result_aggregator reggae_result_aggregator.cpp)
target_link_libraries(reggae_result_aggregator applicationCommon gtcore gtregulation)
set_target_properties(reggae_result_aggregator PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(reggae_result_aggregator)

add_executable(regulator_ora regulator_ora.cpp)
target_link_libraries(regulator_ora applicationCommon gtcore gtregulation)
set_target_properties(regulator_ora PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(regulator_ora)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS reggae reggae_result_aggregator regulator_ora
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
