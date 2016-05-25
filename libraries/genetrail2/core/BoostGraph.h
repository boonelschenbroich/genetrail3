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
#ifndef GT2_CORE_GRAPH_H
#define GT2_CORE_GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <string>

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

typedef boost::adjacency_list< boost::listS, boost::listS, boost::bidirectionalS, boost::property< vertex_identifier_t, std::string >, boost::property< edge_regulation_type_t, std::string > > GraphType;
typedef boost::graph_traits < GraphType >::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits < GraphType >::edge_descriptor edge_descriptor;

#endif //GT2_CORE_GRAPH_H
