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

#ifndef ONTOLOGY_H
#define ONTOLOGY_H


#include <string>
#include <iostream>
#include <boost/date_time/gregorian/gregorian.hpp>

namespace GeneTrail{

	class Ontology{
		
		public:
			
			typedef boost::gregorian::date date_type;
			
			/**
			 * Empty constructor
			 **/
			Ontology() = default;
			
			/**
			 * Default constructor
			 **/
			Ontology(std::string& organism, std::string& url, date_type& creation_date, std::string& clinical_information);
			
			/**
			 * 
			 **/
			Ontology(Ontology&& ontology) = default;
			
			/**
			 * Default copy constructor
			 **/
			Ontology(const Ontology& ontology) = default;
			
			/**
			 * 
			 **/
			std::string organims();
			
			/**
			 * 
			 **/
			void setOrganism(std::string& organims);
			
			/**
			 * 
			 **/
			std::string url();
			
			/**
			 * 
			 **/
			void setUrl(std::string& url);
			
			/**
			 * 
			 **/
			std::string clinicalInformation();
			
			/**
			 * 
			 **/
			void setClinicalInformation(std::string clincial_information);

			/**
			 * 
			 **/
			date_type creationDate();
			
			/**
			 * 
			 **/
			void setCreationDate(boost::gregorian::date creation_date);
			
			/**
			 * 
			 **/
			void writeOntology(std::ostream& out);
			
			/**
			 * 
			 **/
			void readOntology(std::istream& in);

			Ontology& operator=(const Ontology&) = default;
			Ontology& operator=(Ontology&&) = default;
			
		private:
			
			//TODO data_type_ enum????
			
				std::string organism_;
				std::string data_type_;
				std::string url_;
				date_type creation_date_;
				std::string clinical_information_;
			
			
	};
}

#endif // ONTOLOGY_H