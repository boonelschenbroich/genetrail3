#include "BoostGraphProcessor.h"

using namespace GeneTrail;

void BoostGraphProcessor::adeptGraph(GraphType& graph, std::set<std::string> gene_set)
{

	std::set<std::string>::iterator gene_iter;

	typename boost::property_map<GraphType, vertex_identifier_t>::type vertex_ids = boost::get ( vertex_identifier, graph );

	//adapt the graph such that only genes in the testset are cosidered
	boost::graph_traits<GraphType>::vertex_iterator vi, vi_end;
	boost::graph_traits<GraphType>::out_edge_iterator vi3, vi3_end;
	boost::graph_traits<GraphType>::in_edge_iterator vi2, vi2_end;

	int deleted_genes = 0;

	//iterate over all vertices to count the number of nodes in the testset
	for ( tie ( vi, vi_end ) = vertices ( graph ); vi != vi_end; vi++ )
	{
		vertex_descriptor vd = *vi;

		gene_iter = gene_set.find(vertex_ids[vd]);
		bool one_detected = ( gene_iter != gene_set.end() );

		//unify edges, take all in and all outedges of the node as set
		std::set<vertex_descriptor> sources;
		std::set<vertex_descriptor> targets;
		std::set<vertex_descriptor>::iterator sources_iter;
		std::set<vertex_descriptor>::iterator targets_iter;
		std::vector<boost::graph_traits<GraphType>::in_edge_iterator> vit;
		std::vector<boost::graph_traits<GraphType>::out_edge_iterator> vot;

		//if the node has a geneid, add unified edges
		if (!one_detected)
		{
			deleted_genes++;

			for ( tie ( vi2, vi2_end ) = boost::in_edges ( vd, graph ); vi2 != vi2_end; vi2++ )
			{
				sources.insert ( source ( *vi2,graph ) );
				vit.push_back ( vi2 );
			}

			for ( tie ( vi3, vi3_end ) = boost::out_edges ( vd, graph ); vi3 != vi3_end; vi3++ )
			{
				targets.insert ( target ( *vi3,graph ) );
				vot.push_back ( vi3 );
			}

			//remove edges
			for ( unsigned int h = 0; h < vit.size(); h++ )
			{
				boost::remove_edge ( *vit[h], graph );
			}

			for ( unsigned int h = 0; h < vot.size(); h++ )
			{
				vertex_descriptor vd_target = target ( *vot[h], graph );
				vertex_descriptor vd_source = source ( *vot[h], graph );

				if ( vd_source != vd_target )
				{
					boost::remove_edge ( *vot[h], graph );
				}
			}

			//create edges: all sources of inedges of this node to all targets of outedges of this node
			for ( sources_iter = sources.begin(); sources_iter != sources.end(); sources_iter++ )
			{
				for ( targets_iter = targets.begin(); targets_iter != targets.end(); targets_iter++ )
				{
					if ( *sources_iter != *targets_iter && *sources_iter != vd && *targets_iter != vd )
						boost::add_edge ( *sources_iter, *targets_iter, graph );
				}
			}
		}
	}

	for ( tie ( vi, vi_end ) = vertices ( graph ); vi != vi_end; vi++ )
	{
		vertex_descriptor vd = *vi;
		gene_iter = gene_set.find(vertex_ids[vd]);
		bool one_detected = ( gene_iter != gene_set.end() );

		if(!one_detected)
		{
			// This does not work ?!?!!?!!
			//std::cout << "Remove vertex: " << vertex_ids[vd] << std::endl;
			//boost::remove_vertex ( vd, graph );
		}
	}
}

std::set<std::string> BoostGraphProcessor::getVertexSet(GraphType graph)
{
	std::set<std::string> vertex_set;

	typename boost::property_map<GraphType, vertex_identifier_t>::type vertex_ids = boost::get ( vertex_identifier, graph );

	boost::graph_traits<GraphType>::vertex_iterator vi, vi_end;

	for ( tie ( vi, vi_end ) = vertices ( graph ); vi != vi_end; vi++ )
	{
		vertex_descriptor vd = *vi;
		vertex_set.insert(vertex_ids[vd]);
	}

	return vertex_set;
}

