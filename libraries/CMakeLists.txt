project(GENETRAIL2_LIBRARIES_MASTER)

add_subdirectory(genetrail2)

if(GT2_ENABLE_TESTING)
	add_subdirectory(test)
endif()
