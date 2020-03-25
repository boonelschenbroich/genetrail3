#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/EntityDatabase.h>
#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/Exception.h>
#include <genetrail2/core/GeneSet.h>
#include <genetrail2/core/GeneSetReader.h>
#include <genetrail2/core/OverRepresentationAnalysis.h>

#include <genetrail2/enrichment/common.h>
#include <genetrail2/enrichment/EnrichmentAlgorithm.h>
#include <genetrail2/enrichment/SetLevelStatistics.h>
#include <genetrail2/enrichment/Parameters.h>

#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/chrono.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <queue>

using namespace GeneTrail;
using namespace std::chrono_literals;
namespace bpo = boost::program_options;

std::string input, reference, hypothesis, preComputedPValues, method;
bool absolute = false, increasing = false;
int threads = 8;
int already_finished = 0;
int to_do = 0;
boost::mutex mutex_q, mutex_out;
std::queue<std::string> queue;

bool parseArguments(int argc, char* argv[], Params& p)
{
	bpo::variables_map vm;
	bpo::options_description desc;

	addCommonCLIArgs(desc, p);
	desc.add_options()
		("identifiers, d", bpo::value<std::string>(&input)->required(), "A file that stores in each line an identifier file, a tab character, and an output folder. All of identifier files need to contain the same amount of identifiers.")
		("reference, r", bpo::value<std::string>(&reference)->required(), "A file containing identifier line by line.")
		("hypothesis, h", bpo::value<std::string>(&hypothesis)->default_value("two-sided"), "Null hypothesis that should be used.")
		("total-tasks, a", bpo::value<int>(&to_do)->required(), "How many ORAs should be performed.")
		("precomputed-p-values, p", bpo::value<std::string>(&preComputedPValues)->default_value(""), "Precomputed p-values.")
		("threads, t", bpo::value<int>(&threads)->default_value(8), "The number of threads to compute the values. Defaults to 8.")
		("method", bpo::value<std::string>(&method)->default_value("ora"), "The enrichment method to use. Defaults to ora. ");
		
	try
	{
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(),vm);
		bpo::notify(vm);
	}
	catch(bpo::error& e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return false;
	}

	return checkCLIArgs(p);
}

namespace threadNamespace2{
void prepareScores(Scores& scores) {
	if(absolute) {
		std::transform(scores.scores().begin(), scores.scores().end(),
		               scores.scores().begin(),
		               static_cast<double (*)(double)>(std::abs));
	}

	scores.sortByScore(increasing ? Order::Increasing : Order::Decreasing);
}

void thread_job(const DenseMatrix& p_values, const CategoryDBList& cat_list,
				const Category& reference_set, std::shared_ptr<EntityDatabase> db,
				Params p, NullHypothesis& hypothesis_)
{
	while(true){
		std::string line;
		bool empty_queue;
		mutex_q.lock();
		empty_queue = queue.empty();
		if(empty_queue){
			mutex_q.unlock();
			std::this_thread::sleep_for(0.1s);
			mutex_q.lock();
			empty_queue = queue.empty();
		}
		if(!empty_queue){
			line = queue.front();
			queue.pop();
		}
		mutex_q.unlock();
		
		if(empty_queue){
			break;
		}
		
		std::vector<std::string> fields;
		boost::split(fields, line, boost::is_any_of(";"), boost::token_compress_on);
		if(fields.size() < 2) continue;
		if(method == "wilcoxon"){
			p.scores_ = FilePath(fields[0]);
		} else{
			p.identifier_ = FilePath(fields[0]);
		}
		p.out_ = DirectoryPath(fields[1]);
		
		
		GeneSet test_set;
		if(initTestSet(test_set, p) != 0) continue;
		
		try{
			Scores scores(test_set, db);
			if(method == "ora"){
				auto enrichmentAlgorithm = createEnrichmentAlgorithm<Ora>(
					p.pValueMode, reference_set, test_set.toCategory(db, "test"),
					hypothesis_
				);
				run(scores, cat_list, enrichmentAlgorithm, p, true);
			} else if(method == "parallel_ora"){
				auto enrichmentAlgorithm = createEnrichmentAlgorithm<PreprocessedORA>(
					p.pValueMode, reference_set, test_set.toCategory(db, "test"),
					hypothesis_, p_values, p.justScores, p.justPvalues
				);
				run(scores, cat_list, enrichmentAlgorithm, p, true);
			} else if(method == "percentage"){
				auto enrichmentAlgorithm = createEnrichmentAlgorithm<IntersectionPercentage>(
					p.pValueMode, reference_set, test_set.toCategory(db, "test"),
					hypothesis_, p.justScores, p.justPvalues
				);
				run(scores, cat_list, enrichmentAlgorithm, p, true);
			} else if(method == "wilcoxon"){
				auto db = std::make_shared<EntityDatabase>();
				GeneSet test_set;
				CategoryList cat_list;

				if(init(test_set, cat_list, p) != 0) {
					return;
				}
				Scores scores(test_set, db);
				prepareScores(scores);
				
				auto gsea = createEnrichmentAlgorithm<WilcoxonRSTest>(
					p.pValueMode, scores.indices().begin(), scores.indices().end(), increasing ? Order::Increasing : Order::Decreasing);

				run(scores, cat_list, gsea, p, true);
			}
		} catch(IOError& exn){
			std::cerr << "ERROR: Problem loading the output directory " << fields[1];
			std::cerr << " or elements within: " << exn.what() << std::endl;
		}

		boost::lock_guard<boost::mutex> lock{mutex_out};
		std::cout << "Calculating " << ++already_finished << "/" << to_do << " (";
		std::cout << std::fixed << std::setprecision(1) << ((double)already_finished / to_do) * 100.0 << "%)" << std::endl;
	}
	//boost::lock_guard<boost::mutex> lock{mutex_out};
	//std::cout << "I died" << std::endl;
}
}

namespace threadNamespace1{
void input_loader(std::ifstream& input){
	std::string line;
	while(std::getline(input, line)){
		boost::lock_guard<boost::mutex> lock{mutex_q};
		queue.push(line);
	}
}
}

NullHypothesis getHypothesis(const std::string& hypothesis){
	if (hypothesis == "upper-tailed") {
		return NullHypothesis::UPPER_TAILED;
	} else if (hypothesis == "lower-tailed") {
		return NullHypothesis::LOWER_TAILED;
	} else {
		return NullHypothesis::TWO_SIDED;
	}
}

int main(int argc, char* argv[]){
	Params p;
	if(!parseArguments(argc, argv, p)) return -1;
	p.verbose = false;
	
	GeneSet reference_set;
	CategoryList cat_list;
	if(initCategories(cat_list, p) != 0) return -1;

	GeneSetReader reader;
	try{
		reference_set = reader.readGeneList(reference);
	} catch(IOError& exn){
		std::cerr << "ERROR: Failed to read reference set. Reason: " << exn.what() << std::endl;
		return -1;
	}

	auto db = std::make_shared<EntityDatabase>();
	NullHypothesis hypothesis_ = getHypothesis(hypothesis);
	
	DenseMatrix p_values(1,1);
	if(preComputedPValues != ""){
		DenseMatrixReader dm_reader;
		std::ifstream file(preComputedPValues);
		if(!file) {
			std::cerr << "Could not open " << preComputedPValues << " for reading." << std::endl;
			return -1;
		}
		//This is needed since our matrix does not contain any row/col names
		unsigned int opts = DenseMatrixReader::NO_OPTIONS;
		p_values = dm_reader.read(file, opts);
		file.close();
	}
	
	Category ref = reference_set.toCategory(db, "reference");
	
	CategoryDBList category_dbs;
	for(const auto& cat: cat_list){
		try {
			GMTFile gmt(db, cat.second);
			if(!gmt) {
				std::cerr << "WARNING: Could not open database " + cat.first +
									" for reading! Skipping database."
							<< std::endl;
				continue;
			}

			CategoryDatabase b = gmt.read();
			b.setName(cat.first);
			category_dbs.push_back(b);
		} catch(IOError& exn) {
			std::cerr << "WARNING: Could not process category file "
				<< cat.first << "! " << std::endl;
		}
	}
	
	std::ifstream input_strm(input);
	if(!input_strm){
		std::cerr << "Could not open " << input << " for reading." << std::endl;
		return -1;
	}
	std::vector<std::thread> thread_jobs;
	thread_jobs.emplace_back(threadNamespace1::input_loader, boost::ref(input_strm));
	for(int i=0; i < threads; i++){
		thread_jobs.emplace_back(
			threadNamespace2::thread_job,
			boost::cref(p_values), boost::cref(category_dbs),
			boost::cref(ref), db, p, boost::ref(hypothesis_)
		);
	}
	for(auto& t: thread_jobs){
		t.join();
	}
	input_strm.close();
	
	return 0;
}
