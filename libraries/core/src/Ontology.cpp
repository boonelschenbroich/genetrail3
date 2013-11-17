/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#include "Ontology.h"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <vector>
#include <cassert>

 namespace GeneTrail 
{
	
	
	
	Ontology::Ontology(std::string& organism, std::string& url, date_type& creation_date, std::string& clinical_information)
		: organism_(organism),
			url_(url),
			creation_date_(creation_date),
			clinical_information_(clinical_information)
	{
	}
	

	boost::gregorian::date Ontology::creationDate()
	{
		return creation_date_;
	}


	std::string Ontology::clinicalInformation()
	{
		return clinical_information_;
	}

	void Ontology::setClinicalInformation(std::string clincial_information)
	{
		clinical_information_ = clincial_information;
	}

	std::string Ontology::organims()
	{
		return organism_;
	}

	void Ontology::setOrganism(std::string& organims)
	{
		organism_ = organims;
	}

	std::string Ontology::url()
	{
		return url_;
	}

	void Ontology::setUrl(std::string& url)
	{
		url_ = url;
	}

	void Ontology::setCreationDate(boost::gregorian::date creation_date)
	{
		creation_date_ = creation_date;
	}

	void Ontology::readOntology(std::istream& in)
	{
		boost::gregorian::date d;
		std::string clin_inf;
		std::string organism;
		std::string url;
		
		std::string line;
		std::getline(in, line);
		
		boost::trim(line);
		
		std::vector<std::string> splittedLine;
		
		boost::split(splittedLine,line,boost::is_any_of("\t"), boost::token_compress_on);
		assert(splittedLine.size() == 4);
		
		clinical_information_ = splittedLine.at(0);
		
		std::stringstream ss(splittedLine.at(1));
		ss >> creation_date_;
		organism_ = splittedLine.at(2);
		url_ = splittedLine.at(3);
	}

	void Ontology::writeOntology(std::ostream& out)
	{
		out << clinical_information_ << "\t" << creation_date_ << "\t" << organism_ << "\t" << url_;
	}

}