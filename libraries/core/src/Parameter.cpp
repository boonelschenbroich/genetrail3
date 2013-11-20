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


#include "Parameter.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <cassert>

namespace GeneTrail 
{


	Parameter::Parameter(parameter_map row_parameter,parameter_map column_parameter)
		: row_parameter_(std::move(row_parameter)),
		  column_parameter_(std::move(column_parameter))
	{
	}
	
	Parameter::value_type Parameter::col(const Parameter::key_type& key)
	{
		return column_parameter_.at(key);
	}
	
	Parameter::value_type Parameter::row(const Parameter::key_type& key)
	{
		return row_parameter_.at(key);
	}
	
	void Parameter::setCol(const Parameter::key_type& key, Parameter::value_type value)
	{
		column_parameter_.at(key) = value;
	}
	
	void Parameter::setRow(const Parameter::key_type& key, Parameter::value_type value)
	{
		row_parameter_.at(key) = value;
	}
	
	void Parameter::writeParameter(std::ostream& out)
	{
		if(column_parameter_.size() > 0)
		{
			//write parameter of columns
			for(const auto& entry : column_parameter_)
			{
				out << entry.first << "\t" << entry.second 
				<< "\t";
			}
		}
		
		out << "\n";
		
		if(row_parameter_.size() > 0)
		{
			//write parameter of rows
			for(const auto& entry : row_parameter_)
			{
				out << entry.first << "\t" << entry.second << "\t"; 
			}
		}
	}

	void Parameter::readParameter(std::istream& in)
	{
		
// 		assert(column_parameter_.empty());
// 		assert(row_parameter_.empty());
		
		//first the column_parameter_
		std::string line;
		std::getline(in,line);
		
		boost::trim(line);
		
		std::vector<std::string> splittedLine;
		
		boost::split(splittedLine, line, boost::is_any_of("\t"), boost::token_compress_on);
		
		
		if(splittedLine.size() < 2)
		{
			std::cout << "INFO: Parameter: no additional information for columns found" << std::endl;
		}else{
		
			if(splittedLine.size() % 2 != 0)
			{
				std::cerr << "ERROR: Parameter: wrong parameter mapping for columns while reading: " << splittedLine.size() << std::endl;
			}
			
			for(unsigned int i = 0; i < splittedLine.size() ; i+=2)
			{
				key_type first;
				value_type second;
				
				first = splittedLine[i];
				second = splittedLine[i+1];
				
				std::pair<key_type,value_type> pair(first,second);
				
				column_parameter_.insert(pair);
			}
		}
		//row_parameter_
		
		std::getline(in,line);
		boost::trim(line);
		splittedLine.clear();
		boost::split(splittedLine,line,boost::is_any_of("\t"), boost::token_compress_on);
		if(splittedLine.size() < 2)
		{
			std::cout << "INFO: Parameter: no additional information for rows found" << std::endl;
			
		}else{
			if(splittedLine.size() % 2 != 0)
			{
				std::cerr << "ERROR: Parameter: found wrong parameter mapping for rows while reading: " << splittedLine.size() <<  std::endl;
			}
			
			for(unsigned int i = 0; i < splittedLine.size(); i+=2)
			{
				
				key_type first;
				value_type second;
				first = splittedLine[i];
				second = splittedLine[i+1];
				std::pair<key_type,value_type> pair(first,second);
				
				row_parameter_.insert(pair);
			}
		}
	}
	
	void Parameter::insertColumn(const Parameter::key_type& key, Parameter::value_type value)
	{
		column_parameter_.insert(std::make_pair(key,value));
	}

	void Parameter::insertRow(const Parameter::key_type& key, Parameter::value_type value)
	{
		row_parameter_.insert(std::make_pair(key,value));
	}

	Parameter::parameter_map Parameter::col_parameter()
	{
		return column_parameter_;
	}

	Parameter::parameter_map Parameter::row_parameter()
	{
		return row_parameter_;
	}


}