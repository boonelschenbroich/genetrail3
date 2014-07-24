#ifndef GT2_CORE_PERMUTATION_TEST_H
#define GT2_CORE_PERMUTATION_TEST_H

#include <algorithm>
#include <vector>
#include <utility>
#include <random>
#include <iostream>

#include "macros.h"

namespace GeneTrail
{
	template<typename value_type>
	struct TestResult{
		TestResult(std::string n, value_type s, size_t ss)
		:name(n),score(s),counter(0), sampleSize(ss),enriched(false)
		{}

		std::string name;
		value_type score;
		size_t counter;
		size_t sampleSize;
		bool enriched;

		double computePValue(size_t permutations)
		{
			return ((double)counter)/((double)permutations);
		}
	};

	template<typename value_type>
	class GT2_EXPORT PermutationTest
	{
		typedef std::pair<std::string, double> PValue;
		typedef std::vector<PValue> PValues;
		typedef std::vector<TestResult<value_type>> TestResults;

		public:

		PermutationTest() = default;

		PermutationTest(TestResults tests, std::vector<value_type> scores, size_t permutations)
		:permutations_(permutations),tests_(tests),scores_(scores)
		{
			pvalues_.reserve(tests.size());
			// Sort vector of TestResults
			// This makes computation faster
			sort(tests_.begin(), tests_.end(),
			     [](const TestResult<value_type>& a, const TestResult<value_type>& b)
			         ->bool { return a.sampleSize < b.sampleSize; });
		}

		void shuffle(){
			std::mt19937 twister{std::random_device{}()};
			std::shuffle(scores_.begin(), scores_.end(), twister);
		}

		template<typename Function>
		PValues computePValue(Function f){
			for(size_t i=0; i< permutations_; ++i){
				std::cout << "INFO: Running - Permutation test " << i << "/" << permutations_ << std::endl;
				// As we don't want to recompute any values,
				// we save them here.
				shuffle();
				size_t currentSampleSize = -1;
				value_type currentScore = 0.0;
				for(int i=0; i<tests_.size(); ++i){
					// Check if we need to compute new values
					if(tests_[i].sampleSize != currentSampleSize)
					{
						currentSampleSize = tests_[i].sampleSize;
						currentScore = f(scores_.begin(),scores_.begin() + currentSampleSize);
						//std::cout << currentSampleSize << "\t" << currentScore  << std::endl;
					}
					if(tests_[i].enriched)
					{
						if(tests_[i].score <= currentScore){
							++tests_[i].counter;
						}
					}
					else
					{
						if(tests_[i].score >= currentScore){
							++tests_[i].counter;
						}
					}
				}
			}
			for(int i=0; i<tests_.size(); ++i){
				//std::cout << tests_[i].name << "\t" << tests_[i].computePValue(permutations_)  << std::endl;
				pvalues_.push_back(std::make_pair(tests_[i].name, tests_[i].computePValue(permutations_)));
			}
			return pvalues_;
		}

		private:

		size_t permutations_;
		TestResults tests_;
		PValues pvalues_;
		std::vector<value_type> scores_;
	};
}

#endif //GT2_CORE_PERMUTATION_TEST_H

