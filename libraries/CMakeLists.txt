project(GENETRAIL2_LIBRARIES)

if(Boost_FOUND)
	add_subdirectory(core)
else()
	message(STATUS "Boost libraries not found, cannot build any library")
	return()
endif()
