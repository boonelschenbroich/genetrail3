/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 *               2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "FiDePaRunner.h"

using namespace GeneTrail;

void FiDePaRunner::computeDeregulatedPaths(std::string kegg, std::string scores, int pathlength, bool decreasing, bool absolute) const {
    GraphType graph;
    BoostGraphParser graph_parser;
    BoostGraphProcessor graph_processor;
    GeneSetReader reader;

    Pathfinder path_finder;

    graph_parser.readCytoscapeFile<GraphType>(kegg, graph);

    boost::graph_traits<GraphType>::vertex_iterator vid, vid_end;

    auto scoring_file = reader.readNAFile(scores);

	std::vector<GeneSet::Element> sorted_scores;
    if (absolute) {
        sorted_scores = scoring_file.getAbsoluteSortedScores();
    } else {
        sorted_scores = scoring_file.getSortedScores(decreasing);
    }

    std::set<std::string> vertex_set = graph_processor.getVertexSet(graph);

    std::vector<std::string> sorted_gene_list = scoring_file.getIdentifier(scoring_file.intersect(sorted_scores ,vertex_set));

    std::set<std::string> vertex_set_with_scores;

    for (std::string s : sorted_gene_list) {
        vertex_set_with_scores.insert(s);
    }

    graph_processor.adeptGraph(graph, vertex_set_with_scores);

    std::vector<std::vector<std::string> > best_paths;
    std::map<std::string, std::string> regulations;

    pathlength = (pathlength < (signed)sorted_gene_list.size()) ? pathlength : (signed)sorted_gene_list.size();

	std::vector<Path> paths = path_finder.computeDeregulatedPath(graph, sorted_gene_list, pathlength);

	std::vector<std::string> strs;
	boost::split(strs, scores, boost::is_any_of("."));

	int i=2;
	for(auto p : paths){
		p.writeToSIFFile(strs[0] + ".k" + boost::lexical_cast<std::string>(i) + ".sif");
		i = i +1;
	}
}
