#ifndef GRAPH_H
#define GRAPH_H

#ifndef BOOST_GRAPH_UTILITY_HPP
#include <boost/graph/graph_utility.hpp>
#endif

#ifndef BOOST_GRAPH_GRAPH_SEARCH_VISITORS_HPP
#include <boost/graph/visitors.hpp>
#endif

#ifndef BOOST_GRAPH_BREADTH_FIRST_SEARCH_HPP
#include <boost/graph/breadth_first_search.hpp>
#endif

#ifndef BOOST_PROPERTY_MAP_HPP
#include <boost/property_map.hpp>
#endif

#ifndef BOOST_GRAPH_CONNECTED_COMPONENTS_HPP
#include <boost/graph/connected_components.hpp>
#endif

#ifndef BOOST_GRAPH_STRONG_COMPONENTS_HPP
#include <boost/graph/strong_components.hpp>
#endif

#ifndef BOOST_GRAPH_ADJACENCY_LIST_HPP
#include <boost/graph/adjacency_list.hpp>
#endif

#ifndef BOOST_GRAPH_ADJACENCY_MATRIX_HPP
#include <boost/graph/adjacency_matrix.hpp>
#endif

#ifndef BOOST_GRAPH_TRAITS_HPP
#include <boost/graph/graph_traits.hpp>
#endif

#ifndef BOOST_FILTERED_GRAPH_HPP
#include <boost/graph/filtered_graph.hpp>
#endif

#ifndef REVERSE_GRAPH_DWA092300_H_
#include <boost/graph/reverse_graph.hpp>
#endif

#ifndef BOOST_GRAPH_JOHNSON_HPP
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#endif

enum vertex_identifier_t { vertex_identifier };
namespace boost
{
	BOOST_INSTALL_PROPERTY ( vertex, identifier );
}

enum edge_regulation_type_t { edge_regulation_type };
namespace boost
{
	BOOST_INSTALL_PROPERTY ( edge, regulation_type );
}

typedef boost::adjacency_list< boost::listS, boost::vecS, boost::bidirectionalS, boost::property< vertex_identifier_t, std::string >, boost::property< edge_regulation_type_t, std::string > > GraphType;

typedef boost::graph_traits < GraphType >::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits < GraphType >::edge_descriptor edge_descriptor;

#endif
