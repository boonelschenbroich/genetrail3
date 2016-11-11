#include "RegulatorCategoryFileReader.h"

namespace GeneTrail
{
	void RegulatorCategoryFileReader::read_()
	{
		std::ifstream input(file_);
		if(!input) {
			throw GeneTrail::IOError("File (" + file_ +
			                         ") is not open for reading");
		}

		for(std::string line; getline(input, line);) {
			std::vector<std::string> sline(2);
			boost::split(sline, line, boost::is_any_of(" \t"));
			if(sline.size() >= 2) {
				boost::trim(sline[0]);
				boost::trim(sline[1]);
				auto it = categories_.find(sline[0]);
				if(it == categories_.end()) {
					it = categories_.emplace(sline[0], Category(db_)).first;
					it->second.setName(sline[0]);
				}
				it->second.insert(sline[1]);
				reference_.insert(sline[1]);
			} else {
				throw GeneTrail::IOError("Wrong file format.");
			}
		}
	}
}
