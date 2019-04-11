project(Regulation)

####################################################################################################
# Check if the needed libraries are available
####################################################################################################

if(NOT Boost_FOUND)
	message(STATUS "Boost not found, cannot build GeneTrail2 regulation library.")
	return()
endif()

if(NOT RapidJSON_INCLUDE_DIR)
	message(STATUS "RapidJSON_INCLUDE_DIR not specfied, cannot build GeneTrail2 core library.")
	return()
endif()

SET(GT2_CORE_DEP_LIBRARIES
	${Boost_LIBRARIES}
)

if(GENETRAIL2_HAS_GMP)
	list(APPEND GT2_CORE_DEP_LIBRARIES ${GMP_LIBRARIES})
endif()

if(GENETRAIL2_HAS_MPFR)
	list(APPEND GT2_CORE_DEP_LIBRARIES ${MPFR_LIBRARIES})
endif()


####################################################################################################
# Build executable
####################################################################################################

add_library(regulationCommon STATIC common.h common.cpp)
target_include_directories(regulationCommon PRIVATE
	${CMAKE_SOURCE_DIR}/libraries/
	${CMAKE_BINARY_DIR}/libraries/
	"${RapidJSON_INCLUDE_DIR}"
)
GT2_COMPILE_FLAGS(regulationCommon)

add_executable(reggae reggae.cpp)
target_link_libraries(reggae applicationCommon gtcore regulationCommon gtregulation)
set_target_properties(reggae PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(reggae)


add_executable(microReggae microReggae.cpp)
target_link_libraries(microReggae applicationCommon gtcore regulationCommon gtregulation)
set_target_properties(microReggae PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(microReggae)


add_executable(reggae_result_aggregator reggae_result_aggregator.cpp)
target_link_libraries(reggae_result_aggregator applicationCommon gtcore gtregulation)
set_target_properties(reggae_result_aggregator PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(reggae_result_aggregator)

add_executable(regulator_ora regulator_ora.cpp)
target_link_libraries(regulator_ora regulationCommon applicationCommon gtcore gtregulation)
set_target_properties(regulator_ora PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(regulator_ora)

add_executable(regulator_effect_analysis regulator_effect_analysis.cpp)
target_link_libraries(regulator_effect_analysis regulationCommon applicationCommon gtcore gtregulation)
set_target_properties(regulator_effect_analysis PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

configure_file(${CMAKE_SOURCE_DIR}/scripts/regulation/parallel_regulator_ora.in ${CMAKE_SOURCE_DIR}/scripts/regulation/parallel_regulator_ora @ONLY)
file(COPY ${CMAKE_SOURCE_DIR}/scripts/regulation/parallel_regulator_ora
		DESTINATION ${CMAKE_BINARY_DIR}/applications/regulation/
		FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
		GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

GT2_COMPILE_FLAGS(regulator_effect_analysis)

add_executable(csa csa.cpp)
target_link_libraries(csa regulationCommon applicationCommon gtcore gtregulation)
set_target_properties(csa PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

GT2_COMPILE_FLAGS(csa)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS reggae reggae_result_aggregator regulator_ora regulator_effect_analysis csa
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
