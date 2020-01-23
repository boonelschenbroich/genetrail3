/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl <tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_ORA_GROUP_PREFERENCE_H
#define GT2_ORA_GROUP_PREFERENCE_H

#include "macros.h"
#include "DenseMatrix.h"
#include "Metadata.h"
#include <boost/math/distributions/chi_squared.hpp>
#include "multiprecision.h"


#include <utility>
#include <tuple>

#include <boost/math/special_functions/binomial.hpp>

using boost::math::chi_squared;
namespace GeneTrail {

    class GT2_EXPORT ORAGroupPreference {
		public:
			ORAGroupPreference():test_(3){};
			ORAGroupPreference(double threshold): threshold_(threshold), test_(3){};

			/**
			 * This method accepts a (category x sample) matrix file containing p-values
			 * for the category in the respective sample, and a Metadata object with a
			 * key for each sample providing a std::string representing the group of the
			 * corresponding sample. This method calculates a p-value for each category-
			 * group combination on how significant a category is only present in the
			 * respective group.
			 *
			 * @param matrix a (category x sample) matrix with column- and row names.
			 * @param metadata a Metadata object having a key for each sample.
			 * @param result an empty DenseMatrix that is used to store the resulting 
			 *               (category x group) DenseMatrix
			 */
			void calculatePreference(const DenseMatrix& matrix, const Metadata& metadata, DenseMatrix& result) const;
			
			static std::vector<std::string> getGroups(const Metadata& metadata);
			
			static std::vector<std::vector<unsigned int>> getGroupIndices(const DenseMatrix& matrix, const Metadata& metadata, const std::vector<std::string> groups);
			
			static std::vector<std::vector<std::string>> getGroupedSamples(const DenseMatrix& matrix, const Metadata& metadata, const std::vector<std::string> groups);
			
		private:
			unsigned int MINDEX = 0;
			unsigned int LINDEX = 1;
			unsigned int NINDEX = 2;
			unsigned int KINDEX = 3;
			double threshold_;
			chi_squared test_;
			
			void addToTable(std::vector<size_t>& table, double p_value, bool is_current_group) const;
			double computePValue_(const std::vector<size_t>& table) const;
			double computePValue_(size_t m, size_t l, size_t n, size_t k) const;
			
			template <typename T>
			static std::vector<std::vector<T>> getGroupIndices_(
				const DenseMatrix& matrix,
				const Metadata& metadata,
				const std::vector<std::string> groups,
				T dummy
   															)
			{
				std::vector<std::vector<T>> result(groups.size());
				auto col_names = matrix.colNames();
				size_t column_index = -1;
				
				for(const std::string& col_name: col_names){
					column_index++;
					if(metadata.has(col_name)){
						std::string group = get<std::string>(metadata.get(col_name));
						auto elem_it = std::find(groups.begin(), groups.end(), group);
						if(elem_it == groups.end()){
							throw IOError("Internal server error: element " + group + " was " +
							"inserted into the list of groups, but not found later");
						}
						size_t group_index = std::distance(groups.begin(), elem_it);
						result[group_index].push_back(getGroupValue_(column_index, col_name, dummy));
					} else{
						throw IOError("The metadata file had no entry for the sample " + col_name + ".");
					}
				}
				return result;
			}
			
			template <typename T> static T getGroupValue_(size_t, const std::string&, T){return 0;}
			static unsigned int getGroupValue_(size_t c, const std::string&, unsigned int){return c;}
			static std::string getGroupValue_(size_t, const std::string& n, std::string){return n;}
	};
	
}

#endif // GT2_CORE_OVER_REPRESENTATION_ANALYSIS_H

