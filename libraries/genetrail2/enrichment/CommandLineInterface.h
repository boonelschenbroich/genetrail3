#include "PermutationTest.h"

#include <genetrail2/core/macros.h>
#include <genetrail2/core/PValue.h>

#include <boost/optional/optional.hpp>

#include <string>
#include <vector>

namespace boost {
	class any;
	namespace program_options {
		class options_description;
	}
}

namespace GeneTrail {
	class Params;
	class FilePath;
	class DirectoryPath;

	/**
	 * This function adds common arguments to the BOOST commandline parser.
	 *
	 * @param desc BOOST program-options description
	 * @param p Parameter object
	 */
	GT2_EXPORT void addCommonCLIArgs(boost::program_options::options_description& desc, Params& p);

	GT2_EXPORT bool checkCLIArgs(const Params& p);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, PValueMode*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, boost::optional<MultipleTestingCorrection>*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, boost::optional<MatrixHTests>*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, FilePath*,
	              int);

	GT2_EXPORT void validate(boost::any&, const std::vector<std::string>&, DirectoryPath*,
	              int);
}
