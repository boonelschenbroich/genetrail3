/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Tim Kehl <tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORE_BOOST_GRAPH_PARSER_H
#define GT2_CORE_BOOST_GRAPH_PARSER_H

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "BoostGraph.h"
#include "macros.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

namespace GeneTrail {

    class GT2_EXPORT BoostGraphParser {
    public:

        BoostGraphParser() {
        };

        ~BoostGraphParser() {
        };

        /**
         * This method checks if there is a vertex with given identifier contained in the graph.
         * If this vertex is not contained a new vertex is added.
         *
         * @param Identifier string to specify the searched vertex
         * @param Map to find a vertex based on a given identifier
         * @param Boost graph structure
         *
         * @return Pointer to the searched vertex
         */
        template <typename Graph> vertex_descriptor checkForVertexInGraph(std::string identifier, std::map<std::string, vertex_descriptor>& vertex_map, Graph& graph) {
            typename boost::property_map<Graph, vertex_identifier_t>::type vertex_ids = boost::get(vertex_identifier, graph);
            vertex_descriptor vertex_desc;

            if (vertex_map.find(identifier) != vertex_map.end()) {
                vertex_desc = vertex_map.find(identifier)->second;
            } else {
                vertex_desc = boost::add_vertex(graph);
                vertex_ids[vertex_desc] = identifier;
                vertex_map[identifier] = vertex_desc;
            }

            return vertex_desc;
        }

        /**
         * This method constructs a graph from a cytoscape .sif file.
         * FORMAT: ID <tab> REGULATION_TYPE <tab> ID <newline>
         *
         * @param [.sif] file specifying a network structure
         * @param Empty Boost graph structure
         */
        template <typename Graph> void readCytoscapeFile(const std::string& filename, Graph& graph) {
            std::ifstream input_sif;
            std::string current = "";
            std::vector<std::string> entries;

            typename boost::property_map<Graph, edge_regulation_type_t>::type edge_regulations = boost::get(edge_regulation_type, graph);

            vertex_descriptor vertex_desc1, vertex_desc2;
            edge_descriptor edge_desc;

            std::map<std::string, vertex_descriptor> vertex_map;

            //open sif file
            input_sif.open(filename.c_str());

            if (!input_sif) {
                std::cerr << "ERROR: Cannot open: " << filename << std::endl;
            } else {
                while (std::getline(input_sif, current)) {
                    if (current != "") {
                        boost::split(entries, current, boost::is_any_of("\t "));

                        if (entries.size() > 2) {
                            //add vertices
                            boost::trim(entries[0]); //source node
                            boost::trim(entries[1]); //edge regulation
                            boost::trim(entries[2]); //target node

                            vertex_desc1 = checkForVertexInGraph(entries[0], vertex_map, graph);
                            vertex_desc2 = checkForVertexInGraph(entries[2], vertex_map, graph);

                            bool inserted;

                            //add edge
                            boost::tie(edge_desc, inserted) = boost::add_edge(vertex_desc1, vertex_desc2, graph);

                            if (inserted) {
                                edge_regulations[edge_desc] = entries[1];
                            } else {
                                std::cerr << "ERROR: Could not insert edge between " << entries[0] << " and " << entries[2] << filename << std::endl;
                            }
                        }
                    }
                }
            }
        }

        /**
         * This method saves a given graph as a cytoscape .sif file.
         * FORMAT: ID <tab> REGULATION_TYPE <tab> ID <newline>
         *
         * @param Name of the resulting .sif file (without ending)
         * @param Boost graph structure
         */
        template <typename Graph> void writeCytoscapeFile(const std::string& filename, Graph& graph) {
            std::fstream output_sif; //cytoscape sif format
            output_sif.open(filename + ".sif", std::ios::out);

            if (!output_sif) {
                std::cerr << "topology::writeCytoscapeFile -> Cannot open output file: " << filename << " !" << std::endl;
                return;
            }

            typename boost::property_map<Graph, vertex_identifier_t>::type vertex_ids = boost::get(vertex_identifier, graph);
            typename boost::property_map<Graph, edge_regulation_type_t>::type edge_regulations = boost::get(edge_regulation_type, graph);

            typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
            typename boost::graph_traits<Graph>::out_edge_iterator vi2, vi2_end;

            //iterate over all vertices
            for (std::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; vi++) {
                vertex_descriptor vd = *vi;

                //iterate over adjacent vertices of vd
                for (std::tie(vi2, vi2_end) = boost::out_edges(vd, graph); vi2 != vi2_end; vi2++) {
                    output_sif << vertex_ids[vd];

                    if (edge_regulations[*vi2] != "") {
                        output_sif << "\t" << edge_regulations[*vi2] << "\t";
                    } else {
                        output_sif << "\tpp\t";
                    }

                    output_sif << vertex_ids[boost::target(*vi2, graph) ] << std::endl;
                }
            }

            output_sif.close();
        }
    };
}

#endif //GT2_CORE_BOOST_GRAPH_PARSER_H
