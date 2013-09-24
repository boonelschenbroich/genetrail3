project(GENETRAIL2_CORE_LIBRARY_TESTS)

include(../../../test/GeneTrail2TestExecuter.cmake)

####################################################################################################
# Add tests
####################################################################################################

add_gtest(scoring_file_parser_tests LIBRARIES gtcore)
add_gtest(graph_parser_tests        LIBRARIES gtcore)
add_gtest(graph_processor_tests     LIBRARIES gtcore)
add_gtest(pathfinder_tests          LIBRARIES gtcore)
add_gtest(DenseMatrix_tests         LIBRARIES gtcore)
add_gtest(DenseMatrixReader_tests   LIBRARIES gtcore)
add_gtest(DenseMatrixWriter_tests   LIBRARIES gtcore)
