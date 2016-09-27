####################################################################################################
# Path to all files in this directory
####################################################################################################

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "REGULATION")

####################################################################################################
# All files of this directory
####################################################################################################

# We list all files separatly to be able to leave some out, if we do not want to compile them.

# Headers
add_header_to_library(RegulationFileParser.h)
add_header_to_library(RegulationBootstrapper.h)
add_header_to_library(RegulatorGeneAssociationEnrichmentResult.h)
add_header_to_library(RegulatorGeneAssociationEnrichmentAnalysis.h)
add_header_to_library(RegulatorGeneAssociationEnrichmentAlgorithms.h)
add_header_to_library(RegulatorImpactScore.h)

# Sources

add_to_library(RegulatorEnrichmentResultAggregator)