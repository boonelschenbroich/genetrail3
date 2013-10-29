
####################################################################################################
# Path to all files in this directory
####################################################################################################

SET(DIR ${PROJECT_SOURCE_DIR}/src/)

####################################################################################################
# All files of this directory
####################################################################################################

# We list all files separatly to be able to leave some out, if we do not want to compile them.

# Headers
SET(CLUSTER_HEADERS
	${DIR}/METISClusterer.h
	${DIR}/NeighborhoodBuilder.h
	${DIR}/SparseClusterer.h
)

# Sources
SET(CLUSTER_SOURCES
	${DIR}/METISClusterer.cpp
	${DIR}/NeighborhoodBuilder.cpp
	${DIR}/SparseClusterer.cpp
)
