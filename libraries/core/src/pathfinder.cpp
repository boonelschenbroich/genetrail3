#include "pathfinder.h"
#include <algorithm>
#include <map>

using namespace GeneTrail;

// If true the matrix in every step is printed
bool debug=false;

//The number of genes in the list
int numberOfGeneIds;

//Map from gene id to vertex descriptor
typename boost::property_map<GraphType, vertex_identifier_t>::type vertex_ids;

//Saves edge regulations
typename boost::property_map<GraphType, edge_regulation_type_t>::type edge_regs;

//Map from name to rank in gene list
std::map<std::string, int> name2rank;

//Map from rank to vertex_descriptor
std::map<vertex_descriptor,int> vertex_map;

std::vector<vertex_descriptor> nodes;

//Layers of the matrix
std::vector<std::vector<int>> M_1;
std::vector<std::vector<int>> M_2;

//
boost::graph_traits<GraphType>::vertex_iterator vi, vi_end;

// Iterator for all ingoing edges of a vertex
boost::graph_traits<GraphType>::in_edge_iterator iei, iei_end;

// Iterator for all outgoing edges of a vertex
boost::graph_traits<GraphType>::out_edge_iterator oei, oei_end;

// This map contains for each layer a mapping for each vertex to its predecessor
std::vector<std::map<int,int> > best_preds_v;

//Compute RS for all paths of length 1
std::vector<std::vector<int>> running_sums;

void Pathfinder::printMatrix(std::vector<std::vector<int> > m)
{
for(std::vector<int> row : m)
	{
for(int v : row)
		{
			std::string s = (v < 0)?"": " ";
			std::cout << s << v << ", ";
		}

		std::cout << std::endl;
	}
}

int Pathfinder::computeRunningSum(int bpi, int n, int i , int l)
{
	return bpi*n - i*l;
}

void Pathfinder::initializeFields(GraphType& graph, std::vector<std::string>& sorted_gene_list, int& length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string,std::string>& regulations)
{
	//The number of genes in the list
	numberOfGeneIds = sorted_gene_list.size();

	//Map from gene id to vertex descriptor
	vertex_ids = boost::get ( vertex_identifier, graph );
	edge_regs = boost::get ( edge_regulation_type, graph );

	//Map from name to rank in gene list
	for(int i=0; i<numberOfGeneIds; ++i)
	{
		name2rank[sorted_gene_list[i]]=i;
	}

	nodes.resize(numberOfGeneIds);

	for( tie ( vi, vi_end ) = vertices ( graph ); vi != vi_end; vi++)
	{
		vertex_descriptor vd = *vi;
		nodes[name2rank[vertex_ids[vd]]]=vd;
		vertex_map[vd]=name2rank[vertex_ids[vd]];
	}

	// Initialize first layer
	std::vector<std::vector<int>> M( std::vector<std::vector<int> > ( numberOfGeneIds ,std::vector<int> ( numberOfGeneIds,0 ) ) );
	M_1 = M;

	// Fill the first layer
	// This is very simple as the gene list is sorted and we have a bijective mapping
	for(int k=0; k<numberOfGeneIds; ++k)
	{
		for(int i=k; i<numberOfGeneIds; ++i)
		{
			M_1[k][i] = 1;
		}
	}

	// Only needed for debugging
	if(debug)
	{
		printMatrix(M_1);
	}

	// Compute RS for all paths of length 1
	std::vector<std::vector<int>> rs( std::vector<std::vector<int> > ( length ,std::vector<int> ( numberOfGeneIds,0 ) ) );
	running_sums=rs;

	for(int i=0; i<numberOfGeneIds; ++i)
	{
		running_sums[0][i]=computeRunningSum(1,numberOfGeneIds,i+1,1);
	}

	std::cout << "Layer 1" << std::endl;
}

void Pathfinder::findBestPredecessor(GraphType graph, int& best_pred_k, int& best_pred_k_running_sum, int& pred_flag, int l, vertex_descriptor& vd)
{
	for ( ; iei != iei_end; iei++ )
	{
		vertex_descriptor vd_source = source ( *iei, graph );
		int source_kv = vertex_map[vd_source];
		int tmp_rs = running_sums[l-2][source_kv];

		if(pred_flag == 0 || tmp_rs > best_pred_k_running_sum)
		{
			int kv = vertex_map[vd];

			// Check if k is already on the path
			// We have to avoid cycles
			if( (( kv == 0 && M_1[source_kv][kv] == 0 ) || (M_1[source_kv][kv-1] == M_1[source_kv][kv])) && M_1[source_kv][kv] != -1 )
			{
				best_pred_k = source_kv;
				best_pred_k_running_sum = tmp_rs;
				++pred_flag;
			}
		}
	}
}

void Pathfinder::fillNextLayer(int best_pred_k, int k)
{
	for(int i=0; i<numberOfGeneIds; ++i )
	{
		if(k <= i)
		{
			M_2[k][i] = M_1[best_pred_k][i] + 1;
		}
		else
		{
			M_2[k][i] = M_1[best_pred_k][i];
		}
	}
}

void Pathfinder::computeRunningSum(int k, int l, int& max_runnig_sum_k)
{
	int bpi = 1;

	if(M_2[k][0] == 1)
	{
		max_runnig_sum_k = computeRunningSum(bpi,numberOfGeneIds,1,l);
		++bpi;
	}

	for(int i=0; i<numberOfGeneIds-1; ++i)
	{
		if(M_2[k][i] < M_2[k][i+1])
		{
			int running_sum_k = computeRunningSum(bpi,numberOfGeneIds,i+2,l);

			if(bpi == 1 || running_sum_k > max_runnig_sum_k)
			{
				max_runnig_sum_k = running_sum_k;
			}

			++bpi;
		}
	}
}

void Pathfinder::computeDeregulatedPath(GraphType graph, std::vector<std::string> sorted_gene_list, int length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string,std::string>& regulations)
{

	// Initialize fields and first layer of the matrix
	initializeFields(graph, sorted_gene_list, length, best_paths, regulations);

	//Extend the path
	//Fill the layers 2..(length)
	for(int l=2; l <= length; ++l)
	{
		int bestk = -1;
		int bestk_running_sum = -1;

		// Contains a mapping for each vertex to the best predecessor
		std::map<int,int>  best_preds;

		// Initialize second layer
		std::vector<std::vector<int>> M2( std::vector<std::vector<int> > ( numberOfGeneIds ,std::vector<int> ( numberOfGeneIds,-1) ) );
		M_2 = M2;

		for(int k=0; k<numberOfGeneIds; ++k)
		{
			// Find vertex that corresponds to k
			vertex_descriptor vd = nodes[k];

			// Get all ingoing edges
			tie ( iei, iei_end ) = boost::in_edges ( vd, graph );

			// Continue if there are no ingoing edges
			if(iei == iei_end)
			{
				continue;
			}

			// For each target
			// Hold the best values for predecessor
			int best_pred_k;
			int best_pred_k_running_sum;
			int pred_flag=0;

			// Find the best predecessor
			findBestPredecessor(graph, best_pred_k, best_pred_k_running_sum, pred_flag, l , vd);

			// In case there are only cycles possible
			if(pred_flag == 0)
			{
				//std::cout << "only cycles possible" << std::endl;
				continue;
			}

			// Save for each k the best predecessor
			best_preds[k]=best_pred_k;

			// Fill the next layer
			fillNextLayer(best_pred_k,k);


			// Compute the running sum for the current path
			int max_runnig_sum_k;
			computeRunningSum(k, l, max_runnig_sum_k);

			// Save the best running sum
			running_sums[l-1][k]=max_runnig_sum_k;

			if(max_runnig_sum_k > bestk_running_sum)
			{
				bestk = k;
				bestk_running_sum = max_runnig_sum_k;
			}
		}

		// Print all running sums
		if(debug)
		{
			printMatrix(running_sums);
		}

		if(bestk == -1)
		{
			std::cout << "no longer path found" << std::endl;
			break;
		}

		best_preds_v.push_back(best_preds);

		if(debug)
		{
			std::cout << std::endl;
			printMatrix(M_2);
		}

		// Estimate the reversed order
		std::vector<std::string> path;
		std::vector<std::string> rev_path;
		rev_path.push_back(sorted_gene_list[bestk]);

		int tmp_k=bestk;

		for(int i = l-2; i>= 0; --i)
		{
			best_preds = best_preds_v[i];
			tmp_k=best_preds[tmp_k];
			rev_path.push_back(sorted_gene_list[tmp_k]);
		}

		// Revert the path
		// std::vector<std::string> path;
		for(int i=rev_path.size()-1; i >= 0; --i)
		{
			path.push_back(rev_path[i]);

			if(debug)
			{
				std::cout << rev_path[i] << "	" << std::endl;
			}
		}

		// Save all edges with regulations
		for(int i=0; i<(signed)path.size()-1; ++i)
		{
			// get edge regulation
			bool inserted;
			vertex_descriptor vd1 = nodes[name2rank[path[i]]];
			vertex_descriptor vd2 = nodes[name2rank[path[i+1]]];
			edge_descriptor ed;
			boost::tie ( ed, inserted ) = boost::edge ( vd1, vd2, graph );
			regulations[path[i] + path[i+1]]=edge_regs[ed];
		}

		best_paths.push_back(path);

		// As each layer only depends on the layer before we only need to save two layers at once
		M_1=M_2;

		std::cout << "Layer " << l << std::endl;
	}
}
