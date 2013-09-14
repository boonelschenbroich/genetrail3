#ifndef GRAPH_PARSER_H
#define GRAPH_PARSER_H

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "graph.h"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

namespace GeneTrail
{

	class GraphParser
	{
		public:

			GraphParser() {};
			~GraphParser() {};

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
			template <typename Graph> vertex_descriptor checkForVertexInGraph ( std::string identifier, std::map<std::string, vertex_descriptor>& vertex_map, Graph& graph)
			{
				typename boost::property_map<Graph, vertex_identifier_t>::type vertex_ids = boost::get ( vertex_identifier, graph );
				vertex_descriptor vertex_desc;

				if ( vertex_map.find ( identifier ) != vertex_map.end() )
				{
					vertex_desc = vertex_map.find(identifier)->second;
				}
				else
				{
					vertex_desc = boost::add_vertex ( graph );
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
			template <typename Graph> void readCytoscapeFile ( const std::string& filename, Graph& graph)
			{
				std::ifstream input_sif;
				std::string current = "";
				std::vector<std::string> entries;

				typename boost::property_map<Graph, edge_regulation_type_t>::type edge_regulations = boost::get ( edge_regulation_type, graph );

				vertex_descriptor vertex_desc1, vertex_desc2;
				edge_descriptor edge_desc;

				std::map<std::string, vertex_descriptor> vertex_map;

				//open sif file
				input_sif.open ( filename.c_str() );

				if ( !input_sif )
				{
					std::cerr << "readCytoscapeFile -> cannot open: " << filename << std::endl;
				}
				else
				{
					while ( std::getline( input_sif, current ) )
					{
						if ( current != "" )
						{
							boost::split ( entries, current, boost::is_any_of ( "\t " ) );

							if ( entries.size() > 2 )
							{
								//add vertices
								boost::trim ( entries[0] );//source node
								boost::trim ( entries[1] );//edge regulation
								boost::trim ( entries[2] );//target node

								vertex_desc1 = checkForVertexInGraph(entries[0], vertex_map, graph);
								vertex_desc2 = checkForVertexInGraph(entries[2], vertex_map, graph);

								bool inserted;

								//add edge
								boost::tie ( edge_desc, inserted ) = boost::add_edge ( vertex_desc1, vertex_desc2, graph );

								if(inserted)
								{
									edge_regulations[edge_desc]=entries[1];
								}
							}
						}
					}
				}
			}

			///read cytoscape node attribute file and map the node properties in the corresponding property_map
			//format: description<return>
			//node_id1 = label1
			//node_id2 = label2
			//...
			/*template <typename Graph, typename PropertyMap> void mapCytoscapeNodeAttributes ( const std::string& attribute_filename, Graph& graph, PropertyMap& prop )
			{
				std::ifstream input_na;
				std::string current = "";
				std::vector<std::string> entries, labels;
				std::map<std::string, typename PropertyMap::value_type> node2label;
				typename boost::property_map<Graph, vertex_vbnpp_t>::type vindex = boost::get ( vertex_vbnpp, graph );
				typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;

				//open na file
				input_na.open ( attribute_filename.c_str() );

				if ( !input_na )
				{
					std::cerr << "mapCytoscapeNodeAttributes -> cannot open: " << attribute_filename << std::endl;
				}
				else
				{
					while ( pd_.readline ( input_na, current ) )
					{
						if ( current != "" )
						{
							boost::iter_split ( entries, current, boost::first_finder ( " = " ) ); //split line at " = "; first line is description -> is to be ignored

							if ( entries.size() > 1 )
							{
								boost::trim ( entries[0] );
								boost::trim ( entries[1] );
								boost::split ( labels, entries[1], boost::is_any_of ( ";" ) );
								node2label[entries[0]].insert ( node2label[entries[0]].end(), labels.begin(), labels.end() );
							}
						}
					}
				}

				//iterate over vertices and add labels to PropertyMap
				for ( boost::tie ( vi, vi_end ) = boost::vertices ( graph ); vi != vi_end; vi++ )
				{
					vertex_descriptor vd = *vi;
					std::string bnpp_id = vindex[vd];

					if ( node2label.find ( bnpp_id ) != node2label.end() )
					{
						prop[vd] = node2label[bnpp_id];
					}
				}
			}*/

			/**
			 * This method saves a given graph as a cytoscape .sif file.
			 * FORMAT: ID <tab> REGULATION_TYPE <tab> ID <newline>
			 *
			 * @param Name of the resulting .sif file (without ending)
			 * @param Boost graph structure
			 */
			template <typename Graph> void writeCytoscapeFile ( const std::string& filename, Graph& graph )
			{
				std::fstream output_sif; //cytoscape sif format
				output_sif.open ( filename + ".sif", std::ios::out );

				if ( !output_sif)
				{
					std::cerr << "topology::writeCytoscapeFile -> Cannot open output file: " << filename << " !" << std::endl;
					return;
				}

				typename boost::property_map<Graph, vertex_identifier_t>::type vertex_ids = boost::get ( vertex_identifier, graph );
				typename boost::property_map<Graph, edge_regulation_type_t>::type edge_regulations = boost::get ( edge_regulation_type, graph );

				typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
				typename boost::graph_traits<Graph>::out_edge_iterator vi2, vi2_end;

				//iterate over all vertices
				for ( std::tie ( vi, vi_end ) = boost::vertices ( graph ); vi != vi_end; vi++ )
				{
					vertex_descriptor vd = *vi;

					//iterate over adjacent vertices of vd
					for ( std::tie ( vi2, vi2_end ) = boost::out_edges ( vd, graph ); vi2 != vi2_end; vi2++ )
					{
						output_sif << vertex_ids[vd];
						if (edge_regulations[*vi2] != "")
						{
							output_sif << "\t" << edge_regulations[*vi2] << "\t";
						}
						else
						{
							output_sif << "\tpp\t";
						}

						output_sif << vertex_ids[boost::target ( *vi2, graph ) ] << std::endl;
					}
				}

				output_sif.close();
			}
	};
}

#endif

