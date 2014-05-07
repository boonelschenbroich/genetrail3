#ifndef GT2_CORE_GENE_SET_H
#define GT2_CORE_GENE_SET_H

#include "macros.h"
#include "Category.h"

#include <algorithm>
#include <utility>
#include <vector>

namespace GeneTrail
{

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
	class GT2_EXPORT GeneSet
	{
		private:
		typedef std::pair<std::string, value_type> Element;
		typedef std::vector<Element> Container;

		public:
		typedef typename Container::iterator iterator;
		typedef typename Container::const_iterator const_iterator;

		GeneSet() = default;
		GeneSet(const GeneSet&) = default;
		GeneSet(GeneSet&&) = default;

		GeneSet& operator=(const GeneSet&) = default;
		GeneSet& operator=(GeneSet&&) = default;

		iterator begin()
		{
			return container_.begin();
		}

		iterator end()
		{
			return container_.end();
		}

		const_iterator begin() const
		{
			return container_.begin();
		}

		const_iterator end() const
		{
			return container_.end();
		}

        std::vector<std::pair<std::string, value_type> > getScores(){
			return container_;
		}

		void insert(std::pair<std::string, value_type> element)
		{
			container_.push_back(element);
		}

		void insert(std::string element_name, value_type element)
		{
			container_.push_back(std::make_pair(element_name,element));
		}


		int size(){
			return container_.size();
		}

		/**
		 * Getter for the scores object.
		 *
		 * @param Boolean flag indication how to sort the scores (true = decreasing).
		 * @return Sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getSortedScores(bool decreasing){
			std::vector<std::pair<std::string, value_type> > sorted_scores(container_);
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
			std::vector<std::pair<std::string, value_type> > sorted_scores(container_);
			std::sort(sorted_scores.begin(), sorted_scores.end(), increasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getDecreasinglySortedScores(){
			std::vector<std::pair<std::string, value_type> > sorted_scores(container_);
			std::sort(sorted_scores.begin(), sorted_scores.end(), decreasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
        std::vector<std::pair<std::string, value_type> > getAbsoluteSortedScores(){
			std::vector<std::pair<std::string, value_type> > sorted_scores(container_);
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
		std::vector<std::pair<std::string, value_type>> intersect(std::vector<std::pair<std::string, value_type> > scores, std::set<std::string> myset)
		{
			std::vector<std::pair<std::string, value_type>> inter;

			for(auto& elem : scores) {
				std::set<std::string>::iterator setIt;
				setIt = myset.find(elem.first);

				if(setIt != myset.end()) {
					inter.push_back(elem);
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
		std::vector<std::pair<std::string, value_type>> intersect(std::set<std::string> set){
			return intersect(container_, set);
		}

		/**
		 * Computes and returns the intersection with the given set of identifier.
		 *
		 * @return Sorted vector of identifier/score pairs.
		 */
		std::vector<std::pair<std::string, value_type>> sortAndIntersect(std::set<std::string> set, bool decreasing){
			return intersect(getSortedScores(decreasing), set);
		}

		/**
		 * Returns the first k identifier of the given vector.
		 *
		 * @return Vector of identifier.
		 */
        std::vector<std::pair<std::string, value_type>> getFirstK(std::vector<std::pair<std::string, value_type> > scores, int k){
			std::vector<std::pair<std::string, value_type>> firstK;

			for(int i=0; i<k; ++i) {
				firstK.push_back(scores[i]);
			}

			return firstK;
		}

		/**
		 * Returns identifier of the given vector.
		 *
		 * @return Vector of identifier
		 */
		std::vector<std::string> getIdentifier(std::vector<std::pair<std::string, value_type> > scores){
			std::vector<std::string> s;

			for(int i=0; i<scores.size(); ++i) {
				s.push_back(scores[i].first);
			}

			return s;
		}


		const Element& operator[](size_t i) const
		{
			return container_[i];
		}

		Element& operator[](size_t i)
		{
			return container_[i];
		}

		Category toCategory(const std::string& name = "") const
		{
			Category res(name);
			for(const auto& it : container_) {
				res.insert(it.first);
			}
			return res;
		}

		private:
		Container container_;
	};
}

#endif // GT2_CORE_GENE_SET_H

