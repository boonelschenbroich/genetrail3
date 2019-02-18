#############################################################
# Path to all files in this directory
#############################################################

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "CORE")

#############################################################
# All files of this directory
#############################################################

# We list all files separately to be able to leave some out,
# if we do not want to compile them.

# Headers
add_header_to_library(BoostGraph.h)
add_header_to_library(BoostGraphParser.h)
add_header_to_library(CategoryDatabaseFile.h)
add_header_to_library(DenseMatrixIterator.h)
add_header_to_library(GeneSetEnrichmentAnalysis.h)
add_header_to_library(Matrix.h)
add_header_to_library(DependentTTest.h)
add_header_to_library(FTest.h)
add_header_to_library(IndependentTTest.h)
add_header_to_library(OneSampleTTest.h)
add_header_to_library(SignalToNoiseRatio.h)
add_header_to_library(macros.h)
add_header_to_library(MatrixIterator.h)
add_header_to_library(OneSampleWilcoxonSignedRankTest.h)
add_header_to_library(WilcoxonMatchedPairsSignedRankTest.h)
add_header_to_library(WilcoxonRankSumTest.h)
add_header_to_library(AndersonDarlingTest.h)
add_header_to_library(NormalityTest.h)
add_header_to_library(WeightedGeneSetEnrichmentAnalysis.h)
add_header_to_library(ConfidenceInterval.h)
add_header_to_library(BinomialTest.h)
add_header_to_library(NameDatabases.h)
add_header_to_library(MatrixNormalization.h)
add_header_to_library(MatrixTransformation.h)
add_header_to_library(Entropy.h)

# Sources
add_to_library(AbstractMatrix)
add_to_library(BoostGraphProcessor)
add_to_library(Category)
add_to_library(CategoryDatabase)
add_to_library(DenseColumnSubset)
add_to_library(DenseMatrix)
add_to_library(DenseMatrixReader)
add_to_library(DenseMatrixWriter)
add_to_library(DenseRowSubset)
add_to_library(EntityDatabase)
add_to_library(Exception)
add_to_library(File)
add_to_library(FiDePaRunner)
add_to_library(GeneSet)
add_to_library(GeneSetFilters)
add_to_library(GeneSetReader)
add_to_library(GeneSetWriter)
add_to_library(GEO)
add_to_library(GEOGDSParser)
add_to_library(GEOGPLParser)
add_to_library(GEOGSEParser)
add_to_library(GMTFile)
add_to_library(JsonCategoryFile)
add_to_library(MatrixHTest)
add_to_library(MatrixWriter)
add_to_library(Metadata)
add_to_library(misc_algorithms)
add_to_library(Path)
add_to_library(Pathfinder)
add_to_library(PValue)
add_to_library(RMAExpressMatrixReader)
add_to_library(Scores)
add_to_library(SparseMatrix)
add_to_library(SparseMatrixReader)
add_to_library(SparseMatrixWriter)
add_to_library(OverRepresentationAnalysis)
add_to_library(TextFile)
