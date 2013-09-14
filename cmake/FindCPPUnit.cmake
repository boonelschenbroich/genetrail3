#FIND cpp unit test library

#
INCLUDE(LibFindMacros)

#find the headers
FIND_PATH(CPPUnit_INCLUDES NAMES cppunit)

#find the library
FIND_LIBRARY(CPPUnit_LIBRARY NAMES cppunit)

#rename because LIBFIND_PROCESS expects 
SET(CPPUnit_PROCESS_LIBS CPPUnit_LIBRARY)
SET(CPPUnit_PROCESS_INCLUDES CPPUnit_INCLUDES)

#sets CPPUnit_LIBRARIES
LIBFIND_PROCESS(CPPUnit)
