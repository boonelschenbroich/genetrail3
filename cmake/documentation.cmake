file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/doc")

find_package(Doxygen)

if (DOXYGEN_FOUND)
	configure_file("${PROJECT_SOURCE_DIR}/doc/Doxyfile.in" "${PROJECT_BINARY_DIR}/doc/Doxyfile")

	add_custom_target(doc
		COMMAND ${CMAKE_COMMAND} -E echo "Creating html documentation"
		COMMAND ${CMAKE_COMMAND} -E remove_directory doc/html
		COMMAND ${CMAKE_COMMAND} -E chdir doc ${DOXYGEN_EXECUTABLE} Doxyfile
		COMMAND ${CMAKE_COMMAND} -E echo "The documentation has been successfully created."
		COMMAND ${CMAKE_COMMAND} -E echo "You can now open 'doc/html/index.html' in a web browser."
		COMMAND ${CMAKE_COMMAND} -E echo ""
		COMMENT "Build the doxygen documentation"
		VERBATIM
	)
else()
	message(STATUS "Doxygen not found. Disabling all documentation targets!")
endif()
