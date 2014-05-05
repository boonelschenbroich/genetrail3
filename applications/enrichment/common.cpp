#include "common.h"

void addCommonCLIArgs(bpo::options_description& desc, Params& p)
{
	desc.add_options()("help,h", "Display this message")
		("significance,t", bpo::value<double>(&p.significance)->default_value(0.01),"The critical value for rejecting the H0 hypothesis.")
		("categories,c", bpo::value<std::string>(&p.categories)->required(), "A .gmt file containing the categories to be tested.")
		("scores,s", bpo::value<std::string>(&p.scores), "A whitespace seperated file containing identifier and scores.")
		("identifier, i", bpo::value<std::string>(&p.identifier), "A file containing identifier line by line.")
		("minimum,n", bpo::value<int>(&p.minimum)->default_value(0),"Minimum number of genes allowed in categories.")
		("maximum,x", bpo::value<int>(&p.maximum)->default_value(1000),"Maximum number of genes allowed in categories.")
		("output,o", bpo::value<std::string>(&p.out), "Output prefix for text files.")
		("adjustment,a", bpo::value<std::string>(&p.adjustment)->default_value("no"),"P-value adjustment method for multiple testing.");
}

CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat)
{
	//TODO: Stub
	return CategoryList();
}

