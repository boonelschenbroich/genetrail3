#ifndef PATHFINDER_H
#define PATHFINDER_H

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <set>

#include "BoostGraph.h"

namespace GeneTrail
{

	class Pathfinder
	{
		public:

			Pathfinder() {};
			~Pathfinder() {};

			/**
			 * Prints the given matrix
			 *
			 * @param Two dimensional vector representing a matrix
			 */
			void printMatrix(std::vector<std::vector<int> > m);

			/**
			 * Computes the RunningSum (Simplified version of the formula from the FiDePa paper)
			 *
			 * @param Number of genes with rank less or equal to i
			 * @param Number of genes in List
			 * @param The current rank
			 * @param The current length of path
			       * @return The computed Running Sum
			 */
			int computeRunningSum(int bpi, int n, int i , int l);

			/**
			 * Initializes all fields and and the first layer of the matrix
			 *
			 * @param KEGG network as Boost graph structure
			 * @param Gene list sorted by rank
			 * @param Maximal length of path
			 * @param Saves for all k the best path (according to the possibilities of FiDePa)
			 * @param Map containing all edge regulations (GENE_ID + GENE_ID) -> REGULATION
			 */
			void initializeFields(GraphType& graph, std::vector<std::string>& sorted_gene_list, int& length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string,std::string>& regulations);

			/**
			 * Finds the predecessor with best running sum
			 *
			 * @param KEGG network as Boost graph object
			 * @param Saves the index of the best predecessor
			 * @param Saves the running sum for the best predecessor
			 * @param Saves the number of predecessors
			 * @param The current layer
			 * @param The current vertex
			 */
			void findBestPredecessor(GraphType graph, int& best_pred_k, int& best_pred_k_running_sum, int& pred_flag, int l, vertex_descriptor& vd);

			/**
			 * Fills the next layer based on the best predecessor
			 *
			 * @param Best predecessor of k
			 * @param Current vertex
			 */
			void fillNextLayer(int best_pred_k, int k);

			/**
			 * Computes the running sum for a path of length l ending in vertex k
			 *
			 * @param Current vertex
			 * @param Length of the path
			 * @param Saves the running sum
			 */
			void computeRunningSum(int k, int l, int& max_runnig_sum_k);

			/**
			 * Computes best possible deregulated paths according to the FiDePa algorithm.
			 *
			 * @param KEGG network as BOOST graph object
			 * @param Gene list sorted by rank
			 * @param Maximal length of path
			 * @param Saves for all k the best path (according to the possibilities of FiDePa)
			 * @param Map containing all edge regulations (GENE_ID + GENE_ID) -> REGULATION
			 */
			void computeDeregulatedPath(GraphType graph, std::vector<std::string> sorted_gene_list, int length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string,std::string>& regulations);
	};
}

#endif
