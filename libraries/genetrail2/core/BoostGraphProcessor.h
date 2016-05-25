/*
 * GeneTrail2 - An efficient library for interpreting genetic data
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
#ifndef GT2_CORE_BOOST_GRAPH_PROCESSOR_H
#define GT2_CORE_BOOST_GRAPH_PROCESSOR_H

#include "macros.h"

#include "BoostGraph.h"

#include <set>
#include <vector>
#include <iostream>

namespace GeneTrail {

    class GT2_EXPORT BoostGraphProcessor {
    public:

        /**
         * This method adepts a given graph.
         * If vertices of the graph are not contained in the gene set they are deleted and their predecessors and successors are connected.
         *
         * @param graph Boost graph structure
         * @param gene_set set of gene that are allowed in the graph
         */
        void adeptGraph(GraphType& graph, const std::set<std::string>& gene_set);

        /**
         * This method returns a set of vertex identifiers that are contained in the graph.
         *
         * @param graph Boost graph
         * @return Set of vertex identifiers that are contained in the graph
         */
        std::set<std::string> getVertexSet(const GraphType& graph);
    };
}

#endif //GT2_CORE_BOOST_GRAPH_PROCESSOR_H
