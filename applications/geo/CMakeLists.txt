project(geo)

####################################################################################################
# Build executable
####################################################################################################

add_library(common STATIC common.cpp common.h)
target_link_libraries(common gtcore ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(common)
set_target_properties(common PROPERTIES
	INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

add_executable(parseGEO parseGEO.cpp)
target_link_libraries(parseGEO gtcore common ${BOOST_LIBRARIES})
GT2_COMPILE_FLAGS(parseGEO)
set_target_properties(parseGEO PROPERTIES
    INCLUDE_DIRS ${Boost_INCLUDE_DIRS}
)

####################################################################################################
# Build executable
####################################################################################################

install(TARGETS parseGEO
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
)
