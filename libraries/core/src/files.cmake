####################################################################################################
# Path to all files in this directory
####################################################################################################

SET(DIR ${PROJECT_SOURCE_DIR}/src/)

####################################################################################################
# All files of this directory
####################################################################################################

# We list all files separatly to be able to leave some out, if we do not want to compile them. 

# Headers
SET(CORE_HEADERS
	${DIR}/commandline_parser.h
	${DIR}/graph.h
	${DIR}/graph_parser.h
	${DIR}/graph_processor.h
	${DIR}/pathfinder.h
	${DIR}/scoring_file_parser.h
)

# Sources
SET(CORE_SOURCES
	${DIR}/commandline_parser.cpp
	${DIR}/graph_processor.cpp
	${DIR}/pathfinder.cpp
	${DIR}/scoring_file_parser.cpp
)

