/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
 * Copyright (C) 2013 Tobias Frisch <tfrisch@bioinf.uni-sb.de>
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

#ifndef PARAMETER_H
#define PARAMETER_H

#include <iostream>
#include <unordered_map>
#include <string>
//#include <boost/any.hpp>
#include <boost/variant.hpp>


namespace GeneTrail{
	class Parameter{
		
		public:
			
			typedef std::string key_type;
			
			typedef boost::variant<std::string,int,double,long> value_type;
			
			typedef std::unordered_map <key_type, value_type> parameter_map;
			
			/**
			 * Empty constructor
			 **/
			Parameter() = default;
			
			/**
			 * Default construictor
			 **/
			Parameter(const parameter_map& row_parameter,const parameter_map& column_parameter);
			
			Parameter(const Parameter& parameter) = default;
			
			Parameter(Parameter&& parameter) = default;
			
			/**
			 * Return the parameter of a row defined by key
			 **/
			value_type row(const key_type& key);
			
			/**
			 * Return the parameter of a column defined by key
			 **/
			value_type col(const key_type& key);
			
			/**
			 * 
			 **/
			void writeParameter(std::ostream& out);
			
			/**
			 * 
			 **/
			void readParameter(std::istream& in);
			
			/**
			 * Set 
			 **/
			void setRow(const key_type& key, value_type value);
			
			void setCol(const key_type& key, value_type value);
			
			void insertRow(const key_type& key, value_type value);
			
			void insertColumn(const key_type& key, value_type value);
			
			parameter_map row_parameter();
			
			parameter_map col_parameter();
			
			
			Parameter& operator=(const Parameter&) = default;
			Parameter& operator=(Parameter&&) = default;
			
			
		private:
			
			parameter_map row_parameter_;
			
			parameter_map column_parameter_;
			
	};
	
}

#endif // PARAMETER_H