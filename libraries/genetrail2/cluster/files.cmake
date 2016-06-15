
####################################################################################################
# Path to all files in this directory
####################################################################################################

SET(DIR ${PROJECT_SOURCE_DIR})
SET(LIBNAME "CLUSTER")

####################################################################################################
# All files of this directory
####################################################################################################

# We list all files separatly to be able to leave some out, if we do not want to compile them.

add_to_library(METISClusterer)
add_to_library(NeighborhoodBuilder)

