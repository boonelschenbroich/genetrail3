macro(add_header_to_library FILENAME)
	list(APPEND ${LIBNAME}_HEADERS "${DIR}/${FILENAME}")
endmacro()

macro(add_source_to_library FILENAME)
	list(APPEND ${LIBNAME}_SOURCES "${DIR}/${FILENAME}")
endmacro()

macro(add_to_library CLASSNAME)
	add_source_to_library("${CLASSNAME}.cpp")
	add_header_to_library("${CLASSNAME}.h")
endmacro()

