#ifndef GT2_CORE_GENESETFILTER_H
#define GT2_CORE_GENESETFILTER_H

#include "macros.h"

#include <cmath>
#include <string>
#include <utility>

namespace GeneTrail
{
	template<typename T> class GeneSet;

	namespace GeneSetFilter
	{
		/**
		 * Fully virtual base class for GeneSetFilters
		 *
		 * These filters can be applied to a GeneSet and will remove genes
		 * that match certain criteria.
		 */
		class GT2_EXPORT GeneSetFilter
		{
			public:
			// Workaround for GeneSet::Element
			using Element = std::pair<std::string, double>;

			/// Virtual destructor
			virtual ~GeneSetFilter() {}

			/**
			 * This function decides whether an element should be filtered (removed).
			 *
			 * @param e The element that should be tested
			 * @return true iff the element should be removed.
			 */
			virtual bool filter(const Element& e) = 0;

			/**
			 * A preprocessing function that is called before the first element
			 * is filtered.
			 *
			 * @param gene_set The gene set to which the filter will be applied.
			 */
			virtual void setup(const GeneSet<double>& gene_set) = 0;
		};

		/**
		 * A filter that lets every element pass that has a value strictly
		 * larger than a given threshold.
		 */
		class GT2_EXPORT LargerFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default constructor
			 *
			 * @param threshold values smaller or equal than this will be
			 * filtered out.
			 */
			explicit LargerFilter(double threshold) : threshold_(threshold) {}

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return e.second <= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>&) {}

			private:
			double threshold_;
		};

		/**
		 * A filter that lets every element pass that has a value strictly
		 * smaller than a given threshold.
		 */
		class GT2_EXPORT SmallerFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default constructor
			 *
			 * @param threshold values larger or equal than this will be
			 * filtered out.
			 */
			explicit SmallerFilter(double threshold) : threshold_(threshold) {}

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return e.second >= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>&) {}

			private:
			double threshold_;
		};

		/**
		 * A filter that lets every element pass that has an absolute value
		 * strictly larger than a given threshold.
		 */
		class GT2_EXPORT AbsLargerFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default constructor
			 *
			 * @param threshold elements with absolute value strictly smaller
			 * than this will be filtered out.
			 */
			explicit AbsLargerFilter(double threshold) : threshold_(threshold) {}

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return std::fabs(e.second) <= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>&) {}

			private:
			double threshold_;
		};

		/**
		 * A filter that lets every element pass that has an absolute value
		 * strictly smaller than a given threshold.
		 */
		class GT2_EXPORT AbsSmallerFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default constructor
			 *
			 * @param threshold elements with absolute value strictly larger
			 * than this will be filtered out.
			 */
			explicit AbsSmallerFilter(double threshold) : threshold_(threshold) {}

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return std::fabs(e.second) >= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>&) {}

			private:
			double threshold_;
		};

		/**
		 * A filter that lets a predefined fraction of the largest values
		 * pass.
		 */
		class GT2_EXPORT UpperQuantileFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default Constructor
			 *
			 * @param quantile a value in [0,1] that indicates the percentage
			 * of values that should be retained.
			 *
			 * @throws std::out_of_range if the quantile is not in [0,1]
			 */
			explicit UpperQuantileFilter(double quantile);

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return e.second <= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>& gene_set);

			private:
			double quantile_;
			double threshold_;
		};

		/**
		 * A filter that lets a predefined fraction of the smallest values
		 * pass.
		 */
		class GT2_EXPORT LowerQuantileFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default Constructor
			 *
			 * @param quantile a value in [0,1] that indicates the percentage
			 * of values that should be retained.
			 *
			 * @throws std::out_of_range if the quantile is not in [0,1]
			 */
			explicit LowerQuantileFilter(double quantile);

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return e.second >= threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>& gene_set);

			private:
			double quantile_;
			double threshold_;
		};

		/**
		 * A filter that lets a predefined fraction of the largest absolute
		 * values pass.
		 */
		class GT2_EXPORT QuantileFilter : public GeneSetFilter
		{
			public:
			/**
			 * Default Constructor
			 *
			 * @param quantile a value in [0,1] that indicates the percentage
			 * of values that should be retained.
			 *
			 * @throws std::out_of_range if the quantile is not in [0,1]
			 */
			explicit QuantileFilter(double quantile);

			/// @copydoc GeneSetFilter::filter
			bool filter(const Element& e)
			{
				return e.second >= lower_threshold_ &&
				       e.second <= upper_threshold_;
			}

			/// @copydoc GeneSetFilter::setup
			void setup(const GeneSet<double>& gene_set);

			private:
			double quantile_;
			double lower_threshold_;
			double upper_threshold_;
		};
	}
}

#endif // GT2_CORE_GENESETFILTER_H

