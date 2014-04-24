find_package(SWIG 2.0.4)

if(NOT ${SWIG_FOUND})
	message(STATUS "Could not find the SWIG binding generator. Building of language bindings has been disabled")
	return()
endif()

include(${SWIG_USE_FILE})

find_package(PythonLibs)

include_directories("${CMAKE_SOURCE_DIR}/libraries")

add_subdirectory(core)
#add_subdirectory(cluster)
