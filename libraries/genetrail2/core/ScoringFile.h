/*
 * GeneTrail2 - An efficent library for interpreting genetic data
 * Copyright (C) 2014 Tim Kehl tkehl@bioinf.uni-sb.de>
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

#ifndef GT2_CORE_SCORING_FILE_H
#define GT2_CORE_SCORING_FILE_H

#include <vector>
#include <set>
#include <iostream>
#include <utility>
#include <algorithm>

#include "macros.h"
#include "Category.h"

namespace GeneTrail {

    template<typename value_type>
	struct increasing_compare {
        bool operator() (std::pair<std::string, value_type> a, std::pair<std::string, value_type> b) {
            return (a.second < b.second);
        }
    };

	template<typename value_type>
    struct decreasing_compare {
        bool operator() (std::pair<std::string, value_type> a, std::pair<std::string, value_type> b) {
            return (a.second > b.second);
        }
    };

	template<typename value_type>
    struct absolute_compare {
        bool operator() (std::pair<std::string, value_type> a, std::pair<std::string, value_type> b) {
            return (std::abs(a.second) > std::abs(b.second));
        }
    };

	template<typename value_type>
    class GT2_EXPORT ScoringFile {

		public:

		/**
		 * Getter for the scores object.
		 *
		 * @return Vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getScores(){
			return scores_;
		}

		/**
		 * Adds a new score to the scores object.
		 *
		 * @param A new score.
		 */
		void add(std::pair<std::string, value_type> score){
			scores_.push_back(score);
		}

		int size(){
			return scores_.size();
		}

		/**
		 * Getter for the scores object.
		 *
		 * @param Boolean flag indication how to sort the scores (true = decreasing).
		 * @return Sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getSortedScores(bool decreasing){
			std::vector<std::pair<std::string, value_type> > sorted_scores(scores_);
			if (decreasing) {
				std::sort(sorted_scores.begin(), sorted_scores.end(), decreasing_compare<value_type>());
			} else {
				std::sort(sorted_scores.begin(), sorted_scores.end(), increasing_compare<value_type>());
			}
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Increasingly sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getIncreasinglySortedScores(){
			std::vector<std::pair<std::string, value_type> > sorted_scores(scores_);
			std::sort(sorted_scores.begin(), sorted_scores.end(), increasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getDecreasinglySortedScores(){
			std::vector<std::pair<std::string, value_type> > sorted_scores(scores_);
			std::sort(sorted_scores.begin(), sorted_scores.end(), decreasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getAbsoluteSortedScores(){
			std::vector<std::pair<std::string, value_type> > sorted_scores(scores_);
			std::sort(sorted_scores.begin(), sorted_scores.end(), absolute_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Computes and returns the intersecting identifier.
		 *
		 * @param scores
		 * @param set
		 * @return Vector of identifier/score pairs.
	     */
		std::vector<std::string> intersect(std::vector<std::pair<std::string, value_type> > scores, std::set<std::string> myset)
		{
			std::vector<std::string> inter;

			for(auto& elem : scores) {
				std::set<std::string>::iterator setIt;
				setIt = myset.find(elem.first);

				if(setIt != myset.end()) {
					inter.push_back(elem.first);
				}
			}

			return inter;
		}

		/**
		 * Computes and returns the intersecting identifier.
		 *
		 * @param set
		 * @return Vector of identifier/score pairs.
		 */
		std::vector<std::string> intersect(std::set<std::string> set){
			return intersect(scores_, set);
		}

		/**
		 * Computes and returns the intersection with the given set of identifier.
		 *
		 * @return Sorted vector of identifier/score pairs.
		 */
		std::vector<std::string> sortAndIntersect(std::set<std::string> set, bool decreasing){
			return intersect(getSortedScores(decreasing), set);
		}

		/**
		 * Returns the first k identifier of the given vector.
		 *
		 * @return Vector of identifier.
		 */
        std::vector<std::string> getFirstK(std::vector<std::pair<std::string, value_type> > scores, int k){
			std::vector<std::string> firstK;

			for(int i=0; i<k; ++i) {
				firstK.push_back(scores[i].first);
			}

			return firstK;
		}

		/**
		 * Convert to Category.
		 *
		 * @return Category containing all identifer of this ScoringFile.
		 */
		Category convert(std::string name){
			Category c(name);

			for(auto ele : scores_){
				c.insert(ele.first);
			}
			return c;
		}

		/**
		 * Returns identifier of the given vector.
		 *
		 * @return Vector of identifier
		 */
		std::vector<std::string> getIdentifier(std::vector<std::pair<std::string, value_type> > scores){
			return getFirstK(scores, scores.size());
		}

    	protected:

        std::vector<std::pair<std::string, value_type> > scores_;
    };
}

#endif //GT2_CORE_SCORING_FILE_H

