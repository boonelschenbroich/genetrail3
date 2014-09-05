#ifndef GT2_CORE_GENE_SET_H
#define GT2_CORE_GENE_SET_H

#include "macros.h"
#include "Category.h"

#include <algorithm>
#include <utility>
#include <vector>
#include <set>
#include <cmath>

namespace GeneTrail
{
	/**
	 * This struct is used as comparator for sorting
	 */
	template <typename value_type> struct increasing_compare
	{
		bool operator()(std::pair<std::string, value_type> a,
		                std::pair<std::string, value_type> b)
		{
			return (a.second < b.second);
		}
	};

	/**
	 * This struct is used as comparator for sorting
	 */
	template <typename value_type> struct decreasing_compare
	{
		bool operator()(std::pair<std::string, value_type> a,
		                std::pair<std::string, value_type> b)
		{
			return (a.second > b.second);
		}
	};

	/**
	 * This struct is used as comparator for sorting
	 */
	template <typename value_type> struct absolute_compare
	{
		bool operator()(std::pair<std::string, value_type> a,
		                std::pair<std::string, value_type> b)
		{
			return (std::abs(a.second) > std::abs(b.second));
		}
	};

	template<typename value_type>
	class GT2_EXPORT GeneSet
	{
		public:
		typedef std::pair<std::string, value_type> Element;

		private:
		/**
	 	 * Typedefs
	 	 */
		typedef std::vector<Element> Container;

		typedef std::vector<double> _v;
		typedef std::vector<double>::iterator _viter;
		typedef std::function<void(_v, _v)> return_type;

		public:
		typedef typename Container::iterator iterator;
		typedef typename Container::const_iterator const_iterator;

		/**
		 * Constructors
		 */
		GeneSet() = default;
		GeneSet(const GeneSet&) = default;
		GeneSet(GeneSet&&) = default;

		/**
		 * Assignment operators
		 */
		GeneSet& operator=(const GeneSet&) = default;
		GeneSet& operator=(GeneSet&&) = default;

		/**
		 * Begin function
		 *
		 * @return iterator
		 */
		iterator begin()
		{
			return container_.begin();
		}

		/**
		 * End function
		 *
		 * @return iterator
		 */
		iterator end()
		{
			return container_.end();
		}

		/**
		 * Begin function
		 *
		 * @return iterator
		 */
		const_iterator begin() const
		{
			return container_.begin();
		}

		/**
		 * End function
		 *
		 * @return iterator
		 */
		const_iterator end() const
		{
			return container_.end();
		}

		/**
		 * Getter for the container.
		 *
		 * @return Container container_
		 */
		Container getScores(){
			return container_;
		}

		/**
		 * Insert function
		 *
		 * @param The element to insert
		 */
		void insert(std::pair<std::string, value_type> element)
		{
			container_.push_back(element);
		}

		/**
		 * Insert function
		 *
		 * @param element_name Name of the new score
		 * @param element The new score
		 */
		void insert(std::string element_name, value_type element)
		{
			container_.push_back(std::make_pair(element_name,element));
		}

		/**
		 * Gaetter for the size of the container
		 *
		 * @retuen Size of the container_
		 */
		int size(){
			return container_.size();
		}

		/**
		 * Getter for the scores object.
		 *
		 * @param Boolean flag indication how to sort the scores (true = decreasing).
		 * @return Sorted vector of identifier/score pairs.
		 */
		Container getSortedScores(bool decreasing)
		{
			Container sorted_scores(container_);
			if(decreasing) {
				std::sort(sorted_scores.begin(), sorted_scores.end(),
				          decreasing_compare<value_type>());
			} else {
				std::sort(sorted_scores.begin(), sorted_scores.end(),
				          increasing_compare<value_type>());
			}
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Increasingly sorted vector of identifier/score pairs.
		 */
		Container getIncreasinglySortedScores()
		{
			Container sorted_scores(container_);
			std::sort(sorted_scores.begin(), sorted_scores.end(),
			          increasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
		Container getDecreasinglySortedScores()
		{
			Container sorted_scores(container_);
			std::sort(sorted_scores.begin(), sorted_scores.end(),
			          decreasing_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
		Container getAbsoluteSortedScores()
		{
			Container sorted_scores(container_);
			std::sort(sorted_scores.begin(), sorted_scores.end(),
			          absolute_compare<value_type>());
			return sorted_scores;
		}

		/**
		 * Computes and returns the intersecting identifier.
		 *
		 * @param scores
		 * @param set
		 * @return Vector of identifier/score pairs.
	     */
		Container
		intersect(std::vector<std::pair<std::string, value_type>> scores,
		          std::set<std::string> myset)
		{
			Container inter;

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
		Container intersect(std::set<std::string> set){
			return intersect(container_, set);
		}

		/**
		 * Computes and returns the intersection with the given set of identifier.
		 *
		 * @return Sorted vector of identifier/score pairs.
		 */
		Container sortAndIntersect(std::set<std::string> set, bool decreasing)
		{
			return intersect(getSortedScores(decreasing), set);
		}

		/**
		 * Returns the first k identifier of the given vector.
		 *
		 * @return Vector of identifier.
		 */
		Container
		getFirstK(std::vector<std::pair<std::string, value_type>> scores, int k)
		{
			std::vector<std::pair<std::string, value_type>> firstK(k);

			for(int i = 0; i < k; ++i) {
				firstK[i] = scores[i];
			}

			return firstK;
		}

		/**
		 * Returns identifier of the given vector.
		 *
		 * @return Vector of identifier
		 */
		std::vector<std::string>
		getIdentifier(std::vector<std::pair<std::string, value_type>> scores)
		{
			std::vector<std::string> s(scores.size());

			for(size_t i = 0; i < scores.size(); ++i) {
				s[i] = scores[i].first;
			}

			return s;
		}

		std::vector<std::string> getIdentifier()
		{
			return getIdentifier(container_);
		}

		std::vector<std::string> getSortedIdentifier(bool decreasing)
		{
			return getIdentifier(getSortedScores(decreasing));
		}

		std::vector<std::string> getDecreasinglySortedIdentifier()
		{
			return getIdentifier(getSortedScores(true));
		}

		std::vector<std::string> getIncreasinglySortedIdentifier()
		{
			return getIdentifier(getSortedScores(false));
		}

		std::vector<std::string> getAbsoluteSortedIdentifier()
		{
			return getIdentifier(getAbsoluteSortedScores());
		}

		const Element& operator[](size_t i) const
		{
			return container_[i];
		}

		Element& operator[](size_t i)
		{
			return container_[i];
		}

		/**
		 * This function converts the container into a Category
		 *
		 * @param name Name of the created category
		 */
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
	/**
	 * This function applies a given function to all values of the GeneSet.
	 *
	 * @param gene_set
	 * @param f
	 */
	template <typename value_type>
	void transform(GeneSet<value_type>& gene_set, std::function<value_type(value_type)> f)
	{
		for(auto& it : gene_set)
		{
			it.second = f(it.second);
		}
	}

	/**
	 * This function applies the std::abs to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void abs(GeneSet<value_type>& gene_set)
	{
		transform(gene_set, {std::abs<value_type>});
	}

	/**
	 * This function applies the std::sqrt to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void sqrt(GeneSet<value_type>& gene_set)
	{
		transform(gene_set, {static_cast<value_type(*)(value_type)>(std::sqrt)});
	}

	/**
	 * This function applies the std::log to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void log(GeneSet<value_type>& gene_set)
	{
		transform(gene_set, {static_cast<value_type(*)(value_type)>(std::log)});
	}

	/**
	 * This function applies the std::log2 to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void log2(GeneSet<value_type>& gene_set)
	{
		transform(gene_set, {static_cast<value_type(*)(value_type)>(std::log2)});
	}

	/**
	 * This function applies the std::log10 to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void log10(GeneSet<value_type>& gene_set)
	{
		transform(gene_set, {static_cast<value_type(*)(value_type)>(std::log10)});
	}

	/**
	 * This function applies the std::pow to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type>
	void pow(GeneSet<value_type>& gene_set, int n)
	{
		std::transform(gene_set.begin(), gene_set.end(), gene_set.begin(),
		               [n](std::pair<std::string, double>& e) {
			e.second = std::pow(e.second, n);
			return e;
		});
	}

	/**
	 * This function applies the std::pow (power = 2) to all values of the
	 *GeneSet.
	 *
	 * @param gene_set
	 */
	template <typename value_type> void pow2(GeneSet<value_type>& gene_set)
	{
		pow(gene_set, 2);
	}
}

#endif // GT2_CORE_GENE_SET_H

