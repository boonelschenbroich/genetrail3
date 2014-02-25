#include "FiDePaRunner.h"

using namespace GeneTrail;

std::string FiDePaRunner::convertInt(int number) const {
    std::stringstream ss; //create a stringstream
    ss << number; //add number to the stream
    return ss.str(); //return a string with the contents of the stream
}

void FiDePaRunner::writeSifFiles(std::vector<std::vector<std::string> > best_paths, std::map<std::string, std::string> regulations, std::string scores) const {
    std::vector<std::string> strs;
    boost::split(strs, scores, boost::is_any_of("."));
    std::cout << strs[0] << std::endl;
    int ki = 2;

    for (std::vector<std::string> path : best_paths) {
        std::ofstream myfile;
        std::string file = strs[0] + ".k" + convertInt(ki) + ".sif";
        myfile.open(file.c_str());
        ++ki;

        for (unsigned int i = 0; i < path.size() - 1; ++i) {
            auto res = regulations.find(path[i] + path[i + 1]);

            if (res != regulations.end() && res->second != "") {
                myfile << path[i] << "\t" << res->second << "\t" << path[i + 1] << std::endl;
            } else {
                myfile << path[i] << "\tpp\t" << path[i + 1] << std::endl;
            }
        }

        myfile.close();
    }
}

void FiDePaRunner::computeDeregulatedPaths(std::string kegg, std::string scores, int pathlength, bool increasing, bool absolute) const {
    GraphType graph;
    BoostGraphParser graph_parser;
    BoostGraphProcessor graph_processor;
    ScoringFileParser scoring_file_parser;

    Pathfinder path_finder;

    graph_parser.readCytoscapeFile<GraphType>(kegg, graph);

    boost::graph_traits<GraphType>::vertex_iterator vid, vid_end;

    scoring_file_parser.readScoringFile(scores);

    if (absolute) {
        scoring_file_parser.sortScoringFileAbsolute();
    } else {
        scoring_file_parser.sortScoringFile(increasing);
    }

    std::set<std::string> vertex_set = graph_processor.getVertexSet(graph);

    std::vector<std::string> sorted_gene_list = scoring_file_parser.getAllInSet(vertex_set);

    std::set<std::string> vertex_set_with_scores;

    for (std::string s : sorted_gene_list) {
        vertex_set_with_scores.insert(s);
    }

    graph_processor.adeptGraph(graph, vertex_set_with_scores);

    std::vector<std::vector<std::string> > best_paths;
    std::map<std::string, std::string> regulations;

    pathlength = (pathlength < (signed)sorted_gene_list.size()) ? pathlength : (signed)sorted_gene_list.size();

    path_finder.computeDeregulatedPath(graph, sorted_gene_list, pathlength, best_paths, regulations);

    writeSifFiles(best_paths, regulations, scores);
}
