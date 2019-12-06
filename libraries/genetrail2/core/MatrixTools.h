/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2013 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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

#ifndef GT2_MY_MATRIX_TOOLS_H
#define GT2_MY_MATRIX_TOOLS_H

#include "macros.h"
#include "DenseMatrix.h"
#include "DenseColumnSubset.h"

#include <ostream>
#include <vector>

namespace GeneTrail{
	class GT2_EXPORT EmptyGroup : public std::exception{
		public:
		EmptyGroup(const std::string& name) noexcept : groupname_(name) {}

		virtual const char* what() const noexcept
		{
			return (std::string("Group \"") + groupname_ +
			        "\" does not contain any datapoint.").c_str();
		}

		private:
		std::string groupname_;
	};
	
	class GT2_EXPORT MatrixTools{
	public:
		MatrixTools() = default;
		
		std::vector<unsigned int> getIndices(const DenseMatrix& matrix, const std::vector<std::string>& colnames,
											 const std::string& groupname);

		std::tuple<DenseColumnSubset, DenseColumnSubset> splitMatrix(DenseMatrix& matrix,
																	 const std::vector<std::string>& reference,
																	 const std::vector<std::string>& test);
	};
}

#endif //GT2_MATRIX_TOOLS_H
