####################################################################################################
# Path to all files in this directory
####################################################################################################

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "CORE")

####################################################################################################
# All files of this directory
####################################################################################################

# We list all files separatly to be able to leave some out, if we do not want to compile them.

# Headers
add_header_to_library(BoostGraph.h)
add_header_to_library(BoostGraphParser.h)
add_header_to_library(CategoryFile.h)
add_header_to_library(GeneSetEnrichmentAnalysis.h)
add_header_to_library(GeneSetReader.h)
add_header_to_library(Matrix.h)
add_header_to_library(Dataset.h)
add_header_to_library(DenseDataset.h)
add_header_to_library(DenseSubset.h)
add_header_to_library(DependentTTest.h)
add_header_to_library(FTest.h)
add_header_to_library(IndependentTTest.h)
add_header_to_library(WilcoxonMannWhitneyTest.h)
add_header_to_library(PValue.h)
add_header_to_library(ScoringFile.h)
add_header_to_library(macros.h)

# Sources
add_to_library(AbstractMatrix)
add_to_library(BoostGraphProcessor)
add_to_library(Category)
add_to_library(CommandLineParser)
add_to_library(DenseMatrix)
add_to_library(DenseMatrixReader)
add_to_library(DenseMatrixSubset)
add_to_library(DenseMatrixWriter)
add_to_library(Exception)
add_to_library(File)
add_to_library(FiDePaRunner)
add_to_library(GeneSetReader)
add_to_library(GMTFile)
add_to_library(MatrixWriter)
add_to_library(Path)
add_to_library(Pathfinder)
add_to_library(RMAExpressMatrixReader)
add_to_library(SparseMatrix)
add_to_library(SparseMatrixReader)
add_to_library(SparseMatrixWriter)
add_to_library(Ontology)
add_to_library(OverRepresentationAnalysis)
add_to_library(Parameter)
add_to_library(DenseDatasetImpl)
add_to_library(DenseSubset)
