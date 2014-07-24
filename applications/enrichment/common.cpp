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
		("adjustment,a", bpo::value<std::string>(&p.adjustment)->default_value("no"),"P-value adjustment method for multiple testing.")
		("adjust_separately,p", bpo::value(&p.runSeparately)->zero_tokens(),"Indicates if databases are adjusted separatly or combined.");
}

CategoryList getCategoryList(const std::string& catfile_list, const std::string& single_cat)
{
	CategoryList categories;

	std::ifstream input(catfile_list);
	if(!input) {
		throw GeneTrail::IOError("File (" + catfile_list + ") is not open for reading");
	}

	for(std::string line; getline(input, line);) {
		std::vector<std::string> sline(2);
		boost::split(sline, line, boost::is_any_of(" \t"));
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

void readTestSet(GeneSet<double>& test_set, const Params& p)
{
	GeneSetReader<double> reader;
	if(p.scores != "" && p.identifier != ""){
		throw GeneTrail::IOError("Too many input files specified.");
	}else if( p.scores != "" ){
		test_set = reader.readScoringFile(p.scores);
	}else if(p.identifier != ""){
		test_set = reader.readGeneList(p.identifier);
	}else{
		throw GeneTrail::IOError("No input file specified.");
	}
}

PValueList resultVector(const Results& results)
{
	PValueList result;
	result.reserve(results.size());
	for(const auto& jt : results) {
		result.push_back(std::make_pair(jt.first, jt.second->pvalue.convert_to<double>()));
	}
	return result;
}

PValueList resultVector(const AllResults& results)
{
	PValueList result;
	result.reserve(results.size());
	for(const auto& it : results) {
		for(const auto& jt : it.second) {
			result.push_back(std::make_pair(it.first + "\t" + jt.first, jt.second->pvalue.convert_to<double>()));
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

void writeFiles(const std::string& output_dir, const AllResults& all_results)
{

	size_t size = 0;
	for(auto& results_it : all_results)
	{
		size += results_it.second.size();
	}

	for(const auto& database : all_results)
	{
		std::ofstream output(output_dir + "." + database.first + ".txt");
		if(!output)
		{
			throw GeneTrail::IOError("No input file specified.");
		}
		for(const auto& ele : database.second)
		{
			output << ele.second->serialize() << std::endl;
		}
		output.close();
	}
}

int init(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p){

	try{
		readTestSet(test_set, p);
	}
	catch(IOError& exn)
	{
		std::cerr << "ERROR: Failed to read test set. Reason: " << exn.what() << std::endl;
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

std::vector<std::string> getSortedIdentifier(GeneSet<double>& test_set, const Params& p, bool absolute, bool increasing){
	if(p.scores != ""){
		if(absolute){
			return test_set.getAbsoluteSortedIdentifier();
		}else{
			if(increasing){
				return test_set.getIncreasinglySortedIdentifier();
			}else{
				return test_set.getDecreasinglySortedIdentifier();
			}
		}
	}else{
		return test_set.getIdentifier();
	}
}

void updatePValues(Results& results, const PValueList& pvalues)
{
	for(const auto& it : pvalues)
	{
		results[it.first]->pvalue = it.second;
	}
}

void updatePValues(AllResults& results, const PValueList& pvalues)
{
	for(unsigned int i=0; i<pvalues.size(); ++i)
	{
		std::vector<std::string> s;
		boost::split(s, pvalues[i].first, boost::is_any_of("\t"));
		results[s[0]][s[1]]->pvalue = pvalues[i].second;
	}
}

AllResults compute(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p)
{
	AllResults name_to_cat_results;
	for(const auto& cat : cat_list) {
		std::cout << "INFO: Processing - " << cat.first << std::endl;
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
			std::cerr << "WARNING: Could not process category file " << cat.first << "! " << std::endl;
		}
	}

	return name_to_cat_results;
}

void adjustCombined(AllResults& all_results, const Params& p)
{
	auto results = resultVector(all_results);
	results = pvalue<double>::adjustPValues(results, p.adjustment);
	updatePValues(all_results, results);
}

void adjustSeparately(AllResults& all_results, const Params& p)
{
	for(auto& results_it : all_results)
	{
		auto results = resultVector(results_it.second);
		results = pvalue<double>::adjustPValues(results, p.adjustment);
		updatePValues(results_it.second, results);
	}
}

void run(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p, bool computePValue)
{
	AllResults name_to_cat_results(compute(test_set, cat_list, p));
	if(computePValue)
	{
		computePValues(name_to_cat_results);
	}

	// Checks how they should be adjusted
	if(p.runSeparately)
	{
		adjustSeparately(name_to_cat_results, p);
	}
	else
	{
		adjustCombined(name_to_cat_results, p);
	}

	writeFiles(p.out, name_to_cat_results);
}

void run(GeneSet<double>& test_set, CategoryList& cat_list, const Params& p)
{
	run(test_set, cat_list, p, false);
}
