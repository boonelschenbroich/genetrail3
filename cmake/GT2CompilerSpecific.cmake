####################################################################################################
# Compile flags
####################################################################################################

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	LIST(APPEND CXX_FLAGS "-fvisibility=hidden")
	LIST(APPEND CXX_FLAGS "-pedantic")
	LIST(APPEND CXX_FLAGS "-Wall")

	SET(GT2_EXPORT                 "__attribute__((visibility (\"default\")))")
	SET(GT2_LOCAL                  "__attribute__((visibility (\"hidden\")))")
	SET(GT2_EXTERN_VARIABLE "extern __attribute__((visibility (\"default\")))")
else()
	message(STATUS "At the moment, GeneTrail2 supports only the GNU Compiler Collection.")
	return()
endif()

## Enable C++11 mode
if    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	LIST(APPEND CXX_FLAGS "-std=c++11")
	LIST(APPEND CXX_FLAGS "-Wno-deprecated-register")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	LIST(APPEND CXX_FLAGS "-std=c++0x")
endif()

## Glue the compile flags together
string(REPLACE ";" " " GT2_COMPILE_FLAGS "${CXX_FLAGS}")

function(GT2_COMPILE_FLAGS target)
	set_target_properties(${target} PROPERTIES COMPILE_FLAGS ${GT2_COMPILE_FLAGS})
endfunction()
