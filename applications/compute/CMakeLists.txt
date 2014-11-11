project(compute)

####################################################################################################
# Build executable
####################################################################################################

add_executable(computeScores computeScores.cpp)
target_link_libraries(computeScores gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(computeScores)
set_target_properties(computeScores PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(transformScores transformScores.cpp)
target_link_libraries(transformScores gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(transformScores)
set_target_properties(transformScores PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(normalityCheck normalityCheck.cpp)
target_link_libraries(normalityCheck gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(normalityCheck)
set_target_properties(normalityCheck PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Install
####################################################################################################

install(TARGETS computeScores transformScores normalityCheck
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
