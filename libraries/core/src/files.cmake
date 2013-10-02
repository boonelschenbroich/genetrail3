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
	${DIR}/DenseMatrix.h
	${DIR}/DenseMatrixReader.h
	${DIR}/DenseMatrixWriter.h
	${DIR}/CommandLineParser.h
	${DIR}/BoostGraph.h
	${DIR}/BoostGraphParser.h
	${DIR}/BoostGraphProcessor.h
	${DIR}/Pathfinder.h
	${DIR}/ScoringFileParser.h
	${DIR}/RMAExpressMatrixReader.h
)

# Sources
SET(CORE_SOURCES
	${DIR}/DenseMatrix.cpp
	${DIR}/DenseMatrixReader.cpp
	${DIR}/DenseMatrixWriter.cpp
	${DIR}/CommandLineParser.cpp
	${DIR}/BoostGraphProcessor.cpp
	${DIR}/Pathfinder.cpp
	${DIR}/ScoringFileParser.cpp
	${DIR}/RMAExpressMatrixReader.cpp
)

