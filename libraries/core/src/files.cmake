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
	${DIR}/BoostGraph.h
	${DIR}/BoostGraphParser.h
	${DIR}/BoostGraphProcessor.h
	${DIR}/CommandLineParser.h
	${DIR}/DenseMatrix.h
	${DIR}/DenseMatrixReader.h
	${DIR}/DenseMatrixWriter.h
	${DIR}/Matrix.h
	${DIR}/Pathfinder.h
	${DIR}/RMAExpressMatrixReader.h
	${DIR}/ScoringFileParser.h
	${DIR}/SparseMatrix.h
	${DIR}/Ontology.h
	${DIR}/Parameter.h
	${DIR}/Dataset.h
)

# Sources
SET(CORE_SOURCES
	${DIR}/BoostGraphProcessor.cpp
	${DIR}/CommandLineParser.cpp
	${DIR}/DenseMatrix.cpp
	${DIR}/DenseMatrixReader.cpp
	${DIR}/DenseMatrixWriter.cpp
	${DIR}/Matrix.cpp
	${DIR}/Pathfinder.cpp
	${DIR}/RMAExpressMatrixReader.cpp
	${DIR}/ScoringFileParser.cpp
	${DIR}/SparseMatrix.cpp
	${DIR}/Ontology.cpp
	${DIR}/Parameter.cpp
	${DIR}/Dataset.cpp
)

