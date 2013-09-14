project(GENETRAIL2_CORE_LIBRARY_TESTS)

include(../../../test/GeneTrail2TestExecuter.cmake)

####################################################################################################
# Add tests
####################################################################################################

add_gtest(scoring_file_parser_tests SOURCE_FILES ../src/scoring_file_parser.cpp)
add_gtest(graph_parser_tests SOURCE_FILES ../src/graph_parser.h)
add_gtest(graph_processor_tests SOURCE_FILES ../src/graph_processor.cpp ../src/graph_parser.h ../src/scoring_file_parser.cpp)
add_gtest(pathfinder_tests SOURCE_FILES ../src/graph_processor.cpp ../src/graph_parser.h ../src/pathfinder.h)

