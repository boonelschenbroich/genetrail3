project(compute)

####################################################################################################
# Build executable
####################################################################################################


add_executable(computeScores computeScores.cpp)
target_link_libraries(computeScores applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(computeScores)
set_target_properties(computeScores PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(computeNOD nod.cpp)
target_link_libraries(computeNOD applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(computeNOD)
set_target_properties(computeNOD PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(computeNetworkFormular computeNetworkFormular.cpp)
target_link_libraries(computeNetworkFormular applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(computeNetworkFormular)
set_target_properties(computeNetworkFormular PROPERTIES
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

add_executable(normalizeMatrix normalizeMatrix.cpp)
target_link_libraries(normalizeMatrix applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(normalizeMatrix)
set_target_properties(normalizeMatrix PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(transformMatrix transformMatrix.cpp)
target_link_libraries(transformMatrix applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(transformMatrix)
set_target_properties(transformMatrix PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(reorderMatrix reorderMatrix.cpp)
target_link_libraries(reorderMatrix applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(reorderMatrix)
set_target_properties(reorderMatrix PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(NAHandler naHandler.cpp)
target_link_libraries(NAHandler applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(NAHandler)
set_target_properties(NAHandler PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(oraGroupPreference oraGroupPreference.cpp)
target_link_libraries(oraGroupPreference applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(oraGroupPreference)
set_target_properties(oraGroupPreference PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(interactiveUDS interactive_UDS.cpp)
target_link_libraries(interactiveUDS applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(interactiveUDS)
set_target_properties(interactiveUDS PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(computeGroupedScores computeGroupedScores.cpp)
target_link_libraries(computeGroupedScores applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(computeGroupedScores)
set_target_properties(computeGroupedScores PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(precomputedGroupedScores precomputedGroupedScores.cpp)
target_link_libraries(precomputedGroupedScores applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(precomputedGroupedScores)
set_target_properties(precomputedGroupedScores PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(combineReducedEnrichments combineReducedEnrichments.cpp)
target_link_libraries(combineReducedEnrichments applicationCommon gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(combineReducedEnrichments)
set_target_properties(combineReducedEnrichments PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Install
####################################################################################################

install(TARGETS computeScores transformScores normalityCheck normalizeMatrix transformMatrix reorderMatrix interactiveUDS
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
