#include "common.h"

void addCommonCLIArgs(bpo::options_description& desc, Params& p)
{
	desc.add_options()("help,h", "Display this message")
		("significance,t", bpo::value<double>(&p.significance)->default_value(0.01),"The critical value for rejecting the H0 hypothesis.")
		("categories,c", bpo::value<std::string>(&p.categories)->required(), "A .gmt file containing the categories to be tested.")
		("scores,s", bpo::value<std::string>(&p.scores), "A whitespace seperated file containing identifier and scores.")
		("minimum,n", bpo::value<int>(&p.minimum)->default_value(0),"Minimum number of genes allowed in categories.")
		("maximum,x", bpo::value<int>(&p.maximum)->default_value(1000),"Maximum number of genes allowed in categories.")
		("output,o", bpo::value<std::string>(&p.out), "Output prefix for text files.")
		("adjustment,a", bpo::value<std::string>(&p.adjustment)->default_value("no"),"P-value adjustment method for multiple testing.");
}

CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat)
{
	CategoryList categories;
	std::ifstream input(catfile_list);
	if(!input) {
		throw GeneTrail::IOError("File (" + catfile_list + ") is not open for reading");
	}
	for(std::string line; getline(input, line);) {
		std::vector<std::string> sline;
		boost::split(sline, line, boost::is_any_of(" \t"));
		std::cout << line << std::endl;
		if(sline.size() == 2) {
			boost::trim(sline[0]);
			boost::trim(sline[1]);
			categories.emplace_back(std::make_pair(sline[0], sline[1]));
		} else {
			throw GeneTrail::IOError("Wrong file format.");
		}
	}
	return categories;
}


int readTestSet(GeneSet<double>& test_set, const Params& p)
{
	GeneSetReader<double> reader;
	try{
		if( p.scores != "" ){
			test_set = reader.readScoringFile(p.scores);
		}else if(p.identifier != ""){
			test_set = reader.readGeneList(p.identifier);
		}else{
			return -1;
		}
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: Failed to read test set. Reason: " << exn.what() << std::endl;
		return -1;
	}
	return 0;
}


PValueList resultVector(const AllResults& results)
{
	PValueList result;

	for(const auto& it : results) {
		for(const auto& jt : it.second) {
			result.push_back(std::make_pair(it.first + "\t" + jt.first, jt.second.pvalue));
		}
	}

	return result;
}

std::pair<bool, std::pair<int, std::string>> processCategory(Category& c, GeneSet<double>& test_set, const Params& p)
{
	int hits = 0;
	std::string genes = "";
	for(auto s : test_set) {
		if(c.contains(s.first)) {
			genes += genes == "" ? s.first : "," + s.first;
			++hits;
		}
	}
	return std::make_pair(p.minimum <= hits && hits <= p.maximum, std::make_pair(hits, genes));
}

std::map<std::string, std::vector<EnrichmentResult>> splitDatabases(AllResults& all_results, const PValueList& pvalues)
{
	std::map<std::string, std::vector<EnrichmentResult>> result;

	for(const auto& it : pvalues)
	{
		std::vector<std::string> s;
		boost::split(s, it.first, boost::is_any_of("\t"));
		auto find = result.find(s[0]);
		if(find == result.end()){
			result[s[0]] = std::vector<EnrichmentResult>();
		}
		auto ora = all_results[s[0]][s[1]];
		ora.pvalue = it.second;
		result[s[0]].emplace_back(ora);
	}

	return result;
}

void writeFile(const std::string& output_dir, const std::map<std::string,std::vector<EnrichmentResult>>& databases){
	for(const auto& database : databases){
		std::ofstream output;
		output.open (output_dir + "." + database.first + ".txt");
		for(const auto& ele : database.second)
		{
			output << ele.serialize() << std::endl;
		}
		output.close();
	}
}


int init(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p){

	if(readTestSet(test_set, p) != 0){
		return -1;
	}

	if(p.scores != "" && p.identifier != "") {
		std::cerr << "ERROR: Please specify only one input file." << std::endl;
		return -1;
	} else if(p.scores == "" && p.identifier == "") {
		std::cerr << "ERROR: Please specify a input file." << std::endl;
		return -1;
	}

	try
	{
		//TODO: Add single category feature
		cat_list = getCategoryList(p.categories, std::string());
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: Failed to read categories. Reason: " << exn.what() << std::endl;
		return -1;
	}
	return 0;
}

void run(GeneSet<double>& test_set, CategoryList& cat_list, AllResults& name_to_cat_results, const Params& p)
{
	for(const auto& cat : cat_list) {
		std::cout << cat.first << "	" << cat.second << std::endl;
		try
		{
			GMTFile input(cat.second);

			Results name_to_result;
			while(input) {
				Category c = input.read();
				auto pair = processCategory(c, test_set, p);
				if(!pair.first) {
					continue;
				}

				name_to_result.insert(std::make_pair(c.name(), computeEnrichment(c, pair.second)));
			}

			name_to_cat_results.insert(std::make_pair(cat.first, name_to_result));
		}
		catch(IOError& exn)
		{
			std::cerr << "WARNING: Could not process category file " << cat.first << " skipping! " << std::endl;
		}
	}

	auto results = resultVector(name_to_cat_results);
	results = pvalue<double>::adjustPValues(results, p.adjustment);

	writeFile(p.out, splitDatabases(name_to_cat_results, results));
}
