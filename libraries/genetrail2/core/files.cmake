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
add_header_to_library(EnrichmentResult.h)
add_header_to_library(ORAResult.h)
add_header_to_library(GSEAResult.h)
add_header_to_library(CategoryFile.h)
add_header_to_library(GeneSetEnrichmentAnalysis.h)
add_header_to_library(GeneSet.h)
add_header_to_library(GeneSetReader.h)
add_header_to_library(GeneSetWriter.h)
add_header_to_library(Matrix.h)
add_header_to_library(DependentTTest.h)
add_header_to_library(FTest.h)
add_header_to_library(IndependentTTest.h)
add_header_to_library(OneSampleTTest.h)
add_header_to_library(WilcoxonMannWhitneyTest.h)
add_header_to_library(SignalToNoiseRatio.h)
add_header_to_library(PValue.h)
add_header_to_library(macros.h)
add_header_to_library(GeneLabelPermutationTest.h)
add_header_to_library(PermutationTest.h)
add_header_to_library(IndexIterator.h)
add_header_to_library(MatrixHTest.h)
add_header_to_library(OneSampleWilcoxonSignedRankTest.h)
add_header_to_library(WilcoxonMatchedPairsSignedRankTest.h)
add_header_to_library(WilcoxonRankSumTest.h)

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
add_to_library(GEO)
add_to_library(GEOGDSParser)
add_to_library(GEOGPLParser)
add_to_library(GEOGSEParser)
add_to_library(GMTFile)
add_to_library(MatrixWriter)
add_to_library(Path)
add_to_library(Pathfinder)
add_to_library(RMAExpressMatrixReader)
add_to_library(SparseMatrix)
add_to_library(SparseMatrixReader)
add_to_library(SparseMatrixWriter)
add_to_library(OverRepresentationAnalysis)
add_to_library(TextFile)

# Experimental
add_header_to_library(experimental/Dataset.h)
add_header_to_library(experimental/DenseDataset.h)
add_header_to_library(experimental/DenseSubset.h)
add_to_library(experimental/DenseDatasetImpl)
add_to_library(experimental/DenseSubset)
add_to_library(experimental/Ontology)
add_to_library(experimental/Parameter)
