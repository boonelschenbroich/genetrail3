project(GENETRAIL2_REGULATION_LIBRARY_TESTS)

create_test_config_file()

####################################################################################################
# Unit tests for all classes
####################################################################################################

add_gtest(RegulationFile_tests                      LIBRARIES gtcore)
add_gtest(RegulatoryImpactFactors_tests             LIBRARIES gtcore)