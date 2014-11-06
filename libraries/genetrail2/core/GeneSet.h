#ifndef GT2_CORE_GENE_SET_H
#define GT2_CORE_GENE_SET_H

#include "macros.h"

#include <set>
#include <string>
#include <utility>
#include <vector>

namespace GeneTrail
{
	class Category;

	// Forward declaration of GeneSetFilter.
	namespace GeneSetFilter {
		class GeneSetFilter;
	}

	class GT2_EXPORT GeneSet
	{
		/**
	 	 * Typedefs
	 	 */
		public:
		typedef std::pair<std::string, double> Element;

		private:
		typedef std::vector<Element> Container;

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
		iterator begin() { return container_.begin(); }

		/**
		 * End function
		 *
		 * @return iterator
		 */
		iterator end() { return container_.end(); }

		/**
		 * Begin function
		 *
		 * @return iterator
		 */
		const_iterator begin() const { return container_.begin(); }

		/**
		 * End function
		 *
		 * @return iterator
		 */
		const_iterator end() const { return container_.end(); }

		/**
		 * Getter for the container.
		 *
		 * @return Container container_
		 */
		Container getScores() const { return container_; }

		/**
		 * Insert function
		 *
		 * @param The element to insert
		 */
		void insert(const std::pair<std::string, double>& element)
		{
			container_.push_back(element);
		}

		/**
		 * Insert function
		 *
		 * @param element_name Name of the new score
		 * @param element The new score
		 */
		void insert(const std::string& element_name, double element)
		{
			container_.push_back(std::make_pair(element_name,element));
		}

		/**
		 * Getter for the size of the GeneSet
		 *
		 * @return number of elements in the GeneSet
		 */
		int size() const { return container_.size(); }

		/**
		 * Is the GeneSet empty?
		 *
		 * @returns true iff the gene set is empty
		 */
		bool empty() const { return container_.empty(); }

		/**
		 * Getter for the scores object.
		 *
		 * @param Boolean flag indication how to sort the scores (true = decreasing).
		 * @return Sorted vector of identifier/score pairs.
		 */
		Container getSortedScores(bool decreasing) const;

		/**
		 * Getter for the scores object.
		 *
		 * @return Increasingly sorted vector of identifier/score pairs.
		 */
		Container getIncreasinglySortedScores() const;

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
		Container getDecreasinglySortedScores() const;

		/**
		 * Getter for the scores object.
		 *
		 * @return Decreasingly sorted vector of identifier/score pairs.
		 */
		Container getAbsoluteSortedScores() const;

		/**
		 * Computes and returns the intersecting identifier.
		 *
		 * @param scores
		 * @param set
		 * @return Vector of identifier/score pairs.
	     */
		Container intersect(const std::vector<Element>& scores,
		                    const std::set<std::string>& myset) const;

		/**
		 * Computes and returns the intersecting identifier.
		 *
		 * @param set
		 * @return Vector of identifier/score pairs.
		 */
		Container intersect(const std::set<std::string>& set) const;

		/**
		 * Computes and returns the intersection with the given set of identifier.
		 *
		 * @return Sorted vector of identifier/score pairs.
		 */
		Container sortAndIntersect(const std::set<std::string>& set, bool decreasing) const;

		/**
		 * Returns the first k identifier of the given vector.
		 *
		 * @return Vector of identifier.
		 */
		Container
		getFirstK(const std::vector<Element>& scores, int k) const;

		/**
		 * Returns identifier of the given vector.
		 *
		 * @return Vector of identifier
		 */
		std::vector<std::string>
		getIdentifier(const std::vector<Element>& scores) const;

		std::vector<std::string> getIdentifier() const;

		std::vector<std::string> getSortedIdentifier(bool decreasing) const;

		std::vector<std::string> getDecreasinglySortedIdentifier() const;

		std::vector<std::string> getIncreasinglySortedIdentifier() const;

		std::vector<std::string> getAbsoluteSortedIdentifier() const;

		const Element& operator[](size_t i) const { return container_[i]; }

		Element& operator[](size_t i) { return container_[i]; }

		/**
		 * This function converts the container into a Category
		 *
		 * @param name Name of the created category
		 */
		Category toCategory(const std::string& name = "") const;

		/**
		 * Apply a filter to the gene set in order to remove unwanted
		 * values.
		 *
		 * @param gf A pointer to a GeneSetFilter subclass. If nullptr is passed
		 *           nothing happens.
		 *
		 * @return A reference to the gene set. This can be used for chaining.
		 */
		GeneSet& filter(GeneSetFilter::GeneSetFilter* gf);

		/**
		 * Apply a filter to the gene set in order to remove unwanted
		 * values.
		 *
		 * @param gf A GeneSetFilter object.
		 *
		 * @return A reference to the gene set. This can be used for chaining.
		 */
		GeneSet& filter(GeneSetFilter::GeneSetFilter& gf)
		{
			return filter(&gf);
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
	template <typename Func>
	void transform(GeneSet& gene_set, Func f)
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
	void GT2_EXPORT abs(GeneSet& gene_set);

	/**
	 * This function applies the std::sqrt to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT sqrt(GeneSet& gene_set);

	/**
	 * This function applies the std::log to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT log(GeneSet& gene_set);

	/**
	 * This function applies the std::log2 to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT log2(GeneSet& gene_set);

	/**
	 * This function applies the std::log10 to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT log10(GeneSet& gene_set);

	/**
	 * This function applies the std::pow to all values of the GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT pow(GeneSet& gene_set, int n);

	/**
	 * This function applies the std::pow (power = 2) to all values of the
	 *GeneSet.
	 *
	 * @param gene_set
	 */
	void GT2_EXPORT pow2(GeneSet& gene_set);
}

#endif // GT2_CORE_GENE_SET_H

