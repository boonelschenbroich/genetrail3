/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel <dstoeckel@bioinf.uni-sb.de>
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
		 * @param element The element to insert
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
		size_t size() const { return container_.size(); }

		/**
		 * Is the GeneSet empty?
		 *
		 * @returns true iff the gene set is empty
		 */
		bool empty() const { return container_.empty(); }

		/**
		 * Getter for the scores object.
		 *
		 * @param decreasing Boolean flag indication how to sort the scores (true = decreasing).
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
		 * Computes and returns the intersecting identifier between a list
		 * of gene scores and a list of identifiers.
		 *
		 * @param scores A list of identifiers with associated scores.
		 * @param myset A list of identifiers.
		 *
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
		 * Returns the identifier of the given vector of Elements.
		 *
		 * @return A vector of identifier associated with the elements.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string>
		getIdentifier(const std::vector<Element>& scores) const;

		/**
		 * Returns the identifier of GeneSet's elements.
		 *
		 * @return A vector of identifiers contained in the GeneSet.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string> getIdentifier() const;

		/**
		 * Returns a vector of identifiers of the GeneSet's elements.
		 * The identifiers are sorted according to their scores.
		 *
		 * @param decreasing Indicates wheter the identifiers should be sorted
		 *                   increasingly (false) or decreasingly (true).
		 * @return A vector of identifiers contained in the GeneSet.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string> getSortedIdentifier(bool decreasing) const;

		/**
		 * Returns a vector of identifiers of the GeneSet's elements.
		 * The identifiers are sorted decreasingly according to their scores.
		 *
		 * @return A vector of identifiers contained in the GeneSet.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string> getDecreasinglySortedIdentifier() const;

		/**
		 * Returns a vector of identifiers of the GeneSet's elements.
		 * The identifiers are sorted increasingly according to their scores.
		 *
		 * @return A vector of identifiers contained in the GeneSet.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string> getIncreasinglySortedIdentifier() const;

		/**
		 * Returns a vector of identifiers of the GeneSet's elements.
		 * The identifiers are sorted decreasingly according to their
		 * scores' absolute value.
		 *
		 * @return A vector of identifiers contained in the GeneSet.
		 * @warning This function needs to copy the identifers to a new
		 *          vector. Thus this function is very expensive to call!
		 */
		std::vector<std::string> getAbsoluteSortedIdentifier() const;

		/**
		 * Access the i-th element.
		 *
		 * @returns A const reference to the i-th element of the
		 *          GeneSet.
		 */
		const Element& operator[](size_t i) const { return container_[i]; }

		/**
		 * Access the i-th element. Non-const version.
		 *
		 * @returns A const reference to the i-th element of the
		 *          GeneSet.
		 */
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
	 * @param gene_set The GeneSet to which the function should be applied.
	 * @param f A unary function object transforming a real value into another
	 *          real value.
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
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT abs(GeneSet& gene_set);

	/**
	 * This function applies the std::sqrt to all values of the GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT sqrt(GeneSet& gene_set);

	/**
	 * This function applies the std::log to all values of the GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT log(GeneSet& gene_set);

	/**
	 * This function applies the std::log2 to all values of the GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT log2(GeneSet& gene_set);

	/**
	 * This function applies the std::log10 to all values of the GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT log10(GeneSet& gene_set);

	/**
	 * This function applies the std::pow to all values of the GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 * @param n The exponent that should be used in the power function.
	 */
	void GT2_EXPORT pow(GeneSet& gene_set, int n);

	/**
	 * This function applies the std::pow (power = 2) to all values of the
	 * GeneSet.
	 *
	 * @param gene_set The GeneSet to which the function should be applied.
	 */
	void GT2_EXPORT pow2(GeneSet& gene_set);
}

#endif // GT2_CORE_GENE_SET_H

