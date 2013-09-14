#ifndef GRAPH_PROCESSOR_H
#define GRAPH_PROCESSOR_H

#include "graph.h"

#include <set>
#include <vector>
#include <iostream>

namespace GeneTrail
{

	class GraphProcessor
	{
		public:

			/**
			 * This method adepts a given graph.
			 * If vertices of the graph are not contained in the gene set they are deleted and their predecessors and successors are connected.
			 *
			 * @param Boost graph structure
			 * @param set of gene that are allowed in the graph
			 */
			void adeptGraph(GraphType& graph, std::set<std::string> gene_set);

			/**
			 * This method returns a set of vertex identifiers that are contained in the graph.
			 *
			 * @param Boost graph
             * @return Set of vertex identifiers that are contained in the graph
			 */
			std::set<std::string> getVertexSet(GraphType graph);
	};
}

#endif
