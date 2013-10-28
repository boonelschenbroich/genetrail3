project(Tools)

####################################################################################################
# Build executable
####################################################################################################

include_directories(
	"${CMAKE_SOURCE_DIR}/libraries/core/src/"
)

add_executable(mconvert mconvert.cpp)
target_link_libraries(mconvert gtcore)
GT2_COMPILE_FLAGS(mconvert)

INSTALL(TARGETS mconvert
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
