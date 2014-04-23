####################################################################################################
# Compile flags
####################################################################################################

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.4.0")
		message(FATAL_ERROR "GeneTrail 2 requires a clang version >= 3.4")
	endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
		message(FATAL_ERROR "GeneTrail 2 requires a GCC version >= 4.8.0")
	endif()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	LIST(APPEND CXX_FLAGS "-fvisibility=hidden")
	LIST(APPEND CXX_FLAGS "-pedantic")
	LIST(APPEND CXX_FLAGS "-Wall")

	SET(GT2_EXPORT                 "__attribute__((visibility (\"default\")))")
	SET(GT2_LOCAL                  "__attribute__((visibility (\"hidden\")))")
	SET(GT2_EXTERN_VARIABLE "extern __attribute__((visibility (\"default\")))")
else()
	message(STATUS "At the moment, GeneTrail2 supports only Clang and the GNU Compiler Collection.")
	return()
endif()

## Enable C++11 mode
if    ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	LIST(APPEND CXX_FLAGS "-std=c++11")
	LIST(APPEND CXX_FLAGS "-Wno-deprecated-register")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	LIST(APPEND CXX_FLAGS "-std=c++11")
endif()

function(GT2_COMPILE_FLAGS target)
	target_compile_options(${target} PUBLIC ${CXX_FLAGS})
endfunction()
