# Try to find METIS

# This macro will define the following variables

# METIS_FOUND        - METIS is installed on the system
# METIS_INCLUDE_DIRS - The required include directories for METIS
# METIS_LIBRARIES    - The path to the METIS libraries

include(LibFindMacros)

libfind_pkg_check_modules(METIS_PKGCONF metis)

find_path(METIS_INCLUDE_DIR
	NAMES metis.h
	PATHS ${METIS_PKGCONF_INCLUDE_DIRS}
)

find_library(METIS_LIBRARY
	NAMES metis
	PATHS ${METIS_PKGCONF_LIBRARY_DIRS}
)

set(METIS_PROCESS_INCLUDES METIS_INCLUDE_DIR)
set(METIS_PROCESS_LIBS METIS_LIBRARY)
libfind_process(METIS)
