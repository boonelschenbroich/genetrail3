project(GENETRAIL2_LIBRARIES_MASTER)

# Add the libraries directory
add_subdirectory(genetrail2)

if(GT2_ENABLE_TESTING)
	# If testing is enabled add the test subdirectory
	add_subdirectory(test)
endif()
