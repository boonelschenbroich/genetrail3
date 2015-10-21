project(GENETRAIL2_CORE_LIBRARY_TESTS)

create_test_config_file()

####################################################################################################
# Unit tests for all classes
####################################################################################################

add_gtest(BoostGraphParser_tests            LIBRARIES gtcore)
add_gtest(BoostGraphProcessor_tests         LIBRARIES gtcore)
add_gtest(DenseMatrixReader_tests           LIBRARIES gtcore)
add_gtest(DenseMatrixWriter_tests           LIBRARIES gtcore)
add_gtest(DenseMatrix_tests                 LIBRARIES gtcore)
add_gtest(FiDePaRunner_tests                LIBRARIES gtcore)
add_gtest(FishersExactTest_tests            LIBRARIES gtcore)
add_gtest(GMTFile_tests                     LIBRARIES gtcore)
add_gtest(GeneSetEnrichmentAnalysis_tests   LIBRARIES gtcore)
add_gtest(GeneSetReader_tests               LIBRARIES gtcore)
add_gtest(HTests_test                       LIBRARIES gtcore)
add_gtest(HypergeometricTest_tests          LIBRARIES gtcore)
add_gtest(MatrixHTests_tests                LIBRARIES gtcore)
add_gtest(Matrix_tests                      LIBRARIES gtcore)
add_gtest(Metadata_tests                    LIBRARIES gtcore)
add_gtest(MiscAlgorithms_tests              LIBRARIES gtcore)
add_gtest(OverRepresentationAnalysis_tests  LIBRARIES gtcore)
add_gtest(PValue_tests                      LIBRARIES gtcore)
add_gtest(Scores_test                       LIBRARIES gtcore)
add_gtest(Statistic_test                    LIBRARIES gtcore)