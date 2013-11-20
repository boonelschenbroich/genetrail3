/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#ifndef RMAEXPRESSMATRIXREADER_H
#define RMAEXPRESSMATRIXREADER_H

#include "config.h"

#include "DenseMatrixReader.h"

#include <vector>
#include <istream>

namespace GeneTrail
{
	/**
	 * Reader for the RMAExpress binary matrix format
	 *
	 * This class can read the binary files created by
	 * RMAExpressConsole and convert them into our internal
	 * DenseMatrix format.
	 */
	class GT2_EXPORT RMAExpressMatrixReader : public DenseMatrixReader
	{
		public:
			/**
			 * Read a matrix from the given input stream
			 *
			 * \todo The passed options are currently ignored.
			 *
			 * \see DenseMatrixReader::read
			 */
			virtual DenseMatrix read(std::istream& input, unsigned int opts = READ_ROW_NAMES | READ_COL_NAMES) const override;

		private:
			std::string readString_(std::istream& input) const;
			void skipString_(std::istream& input) const;
			void readStringVector_(std::istream& input, std::vector<std::string>& vec) const;
	};

}

#endif // RMAEXPRESSMATRIXREADER_H
