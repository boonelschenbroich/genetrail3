/*
 * GeneTrail2 - An efficent library for interpreting genetic data
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

#ifndef GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H
#define GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H

#include "macros.h"

#include "Category.h"
#include "EnrichmentResult.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <utility>

using namespace boost::multiprecision;

namespace GeneTrail
{
	template<typename float_type, typename big_int_type>
	class GT2_EXPORT GeneSetEnrichmentAnalysis
	{
		constexpr static bool debug = false;

		public:
			int intersectionSize(const Category& category, const std::vector<std::string>& testSet){
				int n = 0;
				for(auto gene : testSet ){
					if(category.contains(gene)){
						++n;
					}
				}
				return n;
			}

			int absMax(int a, int b){
				return std::abs(a) < std::abs(b) ? b : a;
			}

			/**
			 * This method computes the running sum statistic, based on the given categories.
			 *
			 * @param category Category for which the RSc should be computed.
			 * @param testSet Sorted list of genes, based on which the RSc should be computed.
			 * @return The RSc for the given categories.
			 */
			int computeRunningSum(const Category& category, const std::vector<std::string>& testSet){
				int n = testSet.size();
				int l = intersectionSize(category,testSet);
				int nl = n - l;
				int RSc = 0;
				int rs = 0;

				for(const auto& t : testSet){
					if(category.contains(t)){
						rs += nl;
						RSc = absMax(RSc, rs);
					}else{
						rs -= l;
						RSc = absMax(RSc, rs);
					}
				}
				return RSc;
			}

		    std::pair<int, std::vector<std::pair<int, int>>> computeRunningSumP(const Category& category, const std::vector<std::string>& testSet)
		    {
			    int n = testSet.size();
			    int l = intersectionSize(category, testSet);
			    int nl = n - l;
			    int RSc = 0;
			    int rs = 0;

			    std::vector<std::pair<int, int>> points;
				points.push_back(std::make_pair(0,0));
			    bool inc = false;
				int rs_old = 0;
			    for(int i = 0; i < n; ++i) {
				    if(category.contains(testSet[i])) {
						rs += nl;
					    RSc = absMax(RSc, rs);
					    if(rs > rs_old && !inc) {
						    points.push_back(std::make_pair(i-1, rs_old));
						    inc = true;
					    }
				    } else {
					    rs -= l;
					    RSc = absMax(RSc, rs);
					    if(rs < rs_old && inc) {
						    points.push_back(std::make_pair(i-1, rs_old));
						    inc = false;
					    }
				    }
					rs_old = rs;
			    }
				points.push_back(std::make_pair(testSet.size(),0));
			    return std::make_pair(RSc, points);
		    }

		    /**
			 * This method computes the running sum statistic and the corresponding p-value, based on the given categories.
			 *
			 * @param category Category for which the p-value should be computed.
			 * @param testSet Sorted list of genes, based on which the p-value should be computed.
			 * @return The p-value for the given categories.
			 */
			float_type computeTwoSidedPValue(const Category& category, const std::vector<std::string>& testSet){
				int RSc = computeRunningSum(category, testSet);
				return computeTwoSidedPValue(testSet.size(), intersectionSize(category,testSet), RSc);
			}

			/**
			 * This method computes a two-sided p-value for the given running sum statistic via dynamic programming.
			 *
			 * @param n The number of genes in the test set
			 * @param j The number of genes in the category
			 * @param RSc The running sum statistic
			 * @return p-value
			 */
			float_type computeTwoSidedPValue(const int& n, const int& l, const int& RSc){
				return computePValue_(n, l, std::abs(RSc), [](int v, int RSc) {return std::abs(v) < RSc;});
			}


			/**
			 * This method computes a lower-tailed p-value for the given running sum statistic via dynamic programming.
			 *
			 * @param n The number of genes in the test set
			 * @param j The number of genes in the category
			 * @param RSc The running sum statistic
			 * @return p-value
			 */
			float_type computeLeftPValue(const int& n, const int& l, const int& RSc){
				return computePValue_(n, l, std::abs(RSc), [](int v, int RSc) {return -RSc < v;});
			}

			/**
			 *  This method computes a upper-tailed p-value for the given running sum statistic via dynamic programming.
			 *
			 *  @param n The number of genes in the test set
			 *  @param j The number of genes in the category
			 *  @param RSc The running sum statistic
			 *  @return p-value
			 */
			float_type computeRightPValue(const int& n, const int& l, const int& RSc){
				return computePValue_(n, l, std::abs(RSc), [](int v, int RSc) {return v < RSc;});
			}

			template<typename Comparator>
			float_type computePValue_(const int& n, const int& l, const int& RSc, Comparator comp){
				int nl = n-l;

				std::vector<float_type> M(nl+1,0);
				M[0]=1;

				if(debug){
					std::cout << M[0] << " ";
				}

				//Initialize
				for(int k=1; k <= nl; ++k){
				    if(comp(-k * l, RSc)) {
					    M[k] = 1;
					}
					if(debug){
						std::cout << M[k] << " ";
					}
				}

				if(debug){
					std::cout << std::endl;
				}

				for(int i=1; i <= l; ++i){
					int inl = i*nl;
					M[0] = comp(inl, RSc) ? 1 : 0;
					if(debug){
						std::cout << M[0] << " ";
					}
					for(int k=1; k <= nl; ++k){
						if(comp(inl - k*l, RSc)){
							M[k] += M[k-1];
						} else {
							M[k] = 0;
						}
						if(debug){
							std::cout << M[k] << " ";
						}
					}
					if(debug){
						std::cout << std::endl;
					}
				}

				float_type binom = boost::math::binomial_coefficient<float_type>(n, l);

				if(M[nl] >= binom){
					return 0.0;
				}
				return 1.0 - M[nl] / binom;
			}


			/**
			 * This method computes the running sum statistic and the corresponding p-value, based on the given categories.
			 *
			 * @param category Category for which the p-value should be computed.
			 * @param testSet Sorted list of genes, based on which the p-value should be computed.
			 * @return The p-value for the given categories.
			 */
			float_type computeOneSidedPValue(const Category& category, const std::vector<std::string>& testSet){
				int RSc = computeRunningSum(category, testSet);
				if(RSc < 0){
					return computeLeftPValue(testSet.size(), intersectionSize(category,testSet), RSc);
				}else{
					return computeRightPValue(testSet.size(), intersectionSize(category,testSet), RSc);
				}
			}

			std::pair<float_type, std::vector<std::pair<int,int> > > computePValue(const Category& category, const std::vector<std::string>& testSet){
				std::pair<int, std::vector<std::pair<int,int> > > result = computeRunningSumP(category, testSet);
				int RSc = result.first;
			    float_type p = 1.0;
				if(RSc < 0) {
				 	p = computeLeftPValue(testSet.size(), intersectionSize(category, testSet),RSc);
			    } else {
				    p =  computeRightPValue(testSet.size(), intersectionSize(category, testSet),RSc);
			    }
				return std::make_pair(p, result.second);
		    }

			/**
			 * DEPRECATED - Please use the new implementation
			 * This method computes the running sum statistic and the corresponding p-value, based on the given categories.
			 *
			 * @param category Category for which the p-value should be computed.
			 * @param testSet Sorted list of genes, based on which the p-value should be computed.
			 * @return The p-value for the given categories.
			 */
		    float_type computeOneSidedPValueD(const Category& category, const std::vector<std::string>& testSet)
		    {
			    int RSc = computeRunningSum(category, testSet);
				return computeOneSidedPValueD(testSet.size(), intersectionSize(category,testSet), RSc);
		    }

		    /**
			 * DEPRECATED - Please use the new implementation
			 * This method computes a one-sided p-value for the given running sum statistic.
			 * (This is not the original implementation of the Paper)
			 *
			 * @param m The number of genes in the test set
			 * @param l The number of genes in the category
			 * @param RSc The running sum statistic
			 * @return p-value
			 */
			float_type computeOneSidedPValueD(const int& n, const int& l, const int& RSc){
				std::map<int, float_type> pathways;
				pathways[0] = 1;

				bool enriched = RSc >= 0;
				for(int i=0; i<n; i++){
					iterate(pathways, i, RSc, n, l, enriched);
				}

				auto iter = pathways.find(0);

				if(iter != pathways.end()){
					float_type binom = boost::math::binomial_coefficient<float_type>(n, l);
					return 1.0 - iter->second/binom;
				}
				return 1.0;
			}

			/**
			 * This method is needed to fill the map with the number of possible paths.
			 *
			 * @param newpathway Temporary copy of the pathway map
			 * @param value The number of possible paths at this position.
			 * @param newpos The new position in the thought dynamic programming matrix
			 */
			void fillNewPathway(std::map<int,float_type>& newpathway, float_type value, const int& newpos){
				auto iter1 = newpathway.find(newpos);
				if (iter1 == newpathway.end()){
					newpathway[newpos]=value;
				}else{
					newpathway[newpos]=iter1->second + value;
				}
			}

			/**
			 * This method iterates over the previous layer of the thought dynamic programming matrix and fills the current one.
			 *
			 * @param pathway The previous layer of the thought dynamic programming matrix
			 * @param i The index of the current layer
			 * @param RSc The value of the running sum statistic
			 * @param m The number of genes in the test set
			 * @param l The number of genes in the category
			 * @param enriched Boolean flag indicating if we consider a enriched or depleted 
			 */
			void iterate(std::map<int, float_type>& pathway, const int& i, const int& RSc, const int& n, const int& l, const bool& enriched){
				auto iter = pathway.begin();
				std::map<int, float_type> newpathway;

				while(iter != pathway.end()){
					//where could we be in the next step
					//which balls are available
					int black = (iter->first + i*l) / n;
					int white = i - black;

					//having drawn a black ball
					if(black < l){
						int newpos = iter->first + n - l;
						//only relevant if we are smaller RSc
						if(!enriched || (newpos < RSc)){
							fillNewPathway(newpathway, iter->second, newpos);
						}
					}

					//having drawn a white ball
					if(white < (n-l)){
						int newpos = iter->first - l;
					    if(enriched) {
							// checking whether we can still reach zero
							int black = (newpos + (i + 1) * l) / n;
							if((newpos + (l - black) * (n - l) >= 0) && (black <= l)) {
								fillNewPathway(newpathway, iter->second, newpos);
							}
					    } else {
						    if(newpos > RSc){
								fillNewPathway(newpathway, iter->second, newpos);
							}
						}
					}
					++iter;
				}
				pathway = newpathway;
				newpathway.clear();
			}
	};
}

#endif //GT2_CORE_GENE_SET_ENRICHMENT_ANALYSIS_H

