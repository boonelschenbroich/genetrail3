/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2019 Tim Kehl tkehl@bioinf.uni-sb.de>
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
#ifndef GT2_CORE_ENTROPY_H
#define GT2_CORE_ENTROPY_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include <set>
#include <utility>  

#include "macros.h"

#include <boost/numeric/conversion/cast.hpp>

namespace GeneTrail
{

/**
 * A collection of mathematical operations.
 */
namespace Entropy
{
	/**
	* This method calculates the cummulative entropy
	*
	* @param x Sorted vector.
	*/
	template <typename value_type>
	value_type cummulative_entropy_for_sorted_vector(std::vector<value_type>& x)
	{
		size_t m=x.size();
		value_type entropy = value_type(0);
		for(size_t i=1; i<m; ++i){
			value_type i_m = boost::numeric_cast<value_type>(i)/boost::numeric_cast<value_type>(m);
			entropy += (x[i] - x[i-1]) * i_m * log(i_m);
		}
		return -entropy;
	}
	
	/**
	* This method calculates the cummulative entropy
	*
	* @param begin InputIterator.
	* @param end   InputIterator.
	*/
	template <typename value_type, typename InputIterator>
	value_type cummulative_entropy(InputIterator begin, InputIterator end)
	{
		std::vector<value_type> x(std::distance(begin,end));
		std::copy(begin, end, x.begin());
		std::sort(x.begin(), x.end());
		return cummulative_entropy_for_sorted_vector(x);
	}

	/**
	* This method calculates the conditional cummulative entropy h(X|Z).
	*
	* @param X Vector (increasingly sorted).
	* @param Z Vector of indices (increasingly sorted).
	*
	* @return Conditional cummulative entropy h(X|Z)
	*/
	template <typename value_type>
	value_type conditional_cummulative_entropy_for_sorted_vectors(std::vector<value_type>& X, std::vector<size_t>& Z)
	{
		value_type entropy = value_type(0);
		size_t z = Z.size();
		if(z < 2){
			return entropy;
		}
		for(size_t i=1; i<z; ++i){
			value_type i_z = boost::numeric_cast<value_type>(i)/boost::numeric_cast<value_type>(z);
			entropy += (X[Z[i]] - X[Z[i-1]]) * i_z * log(i_z);
		}
		return -entropy;
	}

	/**
	* This method calculates the conditional cummulative entropy h(X|Z).
	*
	* @param X Vector (increasingly sorted).
	* @param Z Vector of vector of indices (both increasingly sorted) - Bins.
	*
	* @return conditional cummulative entropy h(X|Z)
	*/
	template <typename value_type>
	value_type conditional_cummulative_entropy_for_sorted_vectors(std::vector<value_type>& X, std::vector<std::vector<size_t> >& Z)
	{
		value_type entropy = value_type(0);
		size_t m = X.size();
		for(size_t i=0; i<Z.size(); ++i) {
			size_t z = Z[i].size();
			value_type z_m = boost::numeric_cast<value_type>(z)/boost::numeric_cast<value_type>(m);
			entropy += z_m * conditional_cummulative_entropy_for_sorted_vectors(X, Z[i]);
		}
		return entropy;
	}


	struct SquaredDistance
	{
		template <typename value_type> 
		value_type operator()(std::vector<value_type> a, std::vector<value_type> b)
		{
			value_type dist = value_type(0);
			for(size_t i=0; i<a.size(); ++i) {
				value_type diff = (a[i]-b[i]);
				dist += diff * diff;
			}
			return dist;
		}
	};

	template <typename value_type, typename DistanceMeasure>
	std::vector<std::tuple<size_t, size_t, value_type, value_type>> compute_pairwise_distances(const std::vector<value_type>& X, const std::vector<std::vector<value_type> >& Y, DistanceMeasure dist)
	{
		std::vector<std::tuple<size_t, size_t, value_type, value_type>> distances;
		distances.reserve(Y.size()*(Y.size()-1)/2);
		for(size_t i=0; i<Y.size(); ++i) {
			for(size_t j=i+1; j<Y.size(); ++j) {
 				distances.emplace_back(std::make_tuple(i, j, (X[i] - X[j])*(X[i] - X[j]), dist(Y[i], Y[j])));
			}
		}
		std::sort(distances.begin(), distances.end(), [&](const auto& a, const auto& b){
			if(std::get<3>(a) == std::get<3>(b)) {
				return std::get<2>(a) < std::get<2>(b);
			}
			return std::get<3>(a) < std::get<3>(b);
		});
		return std::move(distances);
	}

	template <typename value_type>
	void reorder_vector_inplace(std::vector<std::tuple<size_t, size_t, std::vector<size_t>, value_type>>& vec, size_t start, size_t position, size_t skip) {
		size_t idx = position;
		if(std::get<3>(vec[position]) < std::get<3>(vec[position - 1])) {
			while(std::get<3>(vec[idx]) < std::get<3>(vec[idx-1]) && idx > start) {
				if (idx != skip) {
					std::swap(vec[idx], vec[idx-1]);
				}
				idx -= 1;
			}
		} else if (std::get<3>(vec[position]) > std::get<3>(vec[position - 1])) {
			while(std::get<3>(vec[idx]) > std::get<3>(vec[idx+1]) && idx < vec.size()-1) {
				if(idx != skip) {
					std::swap(vec[idx], vec[idx+1]);
				}
				idx += 1;
			}
		}
	}

	/**
	* This method calculates the conditional cummulative entropy h(X|Y).
	*
	* @param X Vector of scores.
	* @param Y Vector of vector of scores.
	*
	* @return conditional cummulative entropy h(X|Y)
	*/
	template <typename value_type>
	std::vector<value_type> conditional_cummulative_entropy_estimator_for_sorted_vectors(std::vector<value_type>& X, std::vector<std::vector<value_type> >& Y)
	{
		//Sort indices according to Y
		SquaredDistance dist;
		std::vector<size_t> indices(X.size());	
		std::vector<value_type> nullv(Y[0].size(), value_type(0));
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b){
			auto da = dist(Y[a], nullv); 
			auto db = dist(Y[b], nullv);
			return da < db;
		});
		std::vector<std::vector<value_type>> Y_sorted(Y.size());
		std::copy(Y.begin(), Y.end(), Y_sorted.begin());
		std::sort(Y_sorted.begin(), Y_sorted.end(), [&](auto a, auto b){
			auto da = dist(a, nullv); 
			auto db = dist(b, nullv);
			return da < db;
		});

		//Initial bins
		std::vector<value_type> bin_entropy(Y.size());
		for(size_t i=0; i<Y_sorted.size(); ++i) {
			bin_entropy[i] = value_type(0);
		}

		//Calculate entropy of neighbours
		std::vector<std::tuple<size_t, size_t, std::vector<size_t>, value_type>> entropy_between_neighbors(Y_sorted.size()-1);
		for(size_t i=1; i<Y_sorted.size(); ++i){
			std::vector<size_t> nbins(2);
			nbins[0] = indices[i];
			nbins[1] = indices[i-1];
			std::sort(nbins.begin(), nbins.end());
			value_type z_n = boost::numeric_cast<value_type>(nbins.size()) / boost::numeric_cast<value_type>(X.size());
			entropy_between_neighbors[i-1] = std::make_tuple(i-1, i, nbins, z_n * conditional_cummulative_entropy_for_sorted_vectors(X, nbins));
		}

		std::sort(entropy_between_neighbors.begin(), entropy_between_neighbors.end(), [&](const auto& a, const auto& b){
			return std::get<3>(a) < std::get<3>(b);
		});

		/*std::cout << std::endl << "##############################################" << std::endl;
		for(auto& e: entropy_between_neighbors) {
			std::cout << "Entropy(" << std::get<0>(e) << ", " << std::get<1>(e) << ") - [";
			for(size_t i=0; i<std::get<2>(e).size(); ++i) {
				std::cout << std::get<2>(e)[i] << ", ";
			}
			std::cout << "] - " << std::get<3>(e) << std::endl;
		}*/

		std::string tmp;
		value_type entropy = 0.0;
		std::vector<value_type> entropies;
		entropies.reserve(X.size());
		entropies.emplace_back(value_type(0));
		for(size_t i=0; i<entropy_between_neighbors.size(); ++i){
			size_t left = std::get<0>(entropy_between_neighbors[i]);
			size_t right = std::get<1>(entropy_between_neighbors[i]);
			
			//Update entropy
			//entropy -= bin_entropy[left];
			//entropy -= bin_entropy[right];
			entropy += std::get<3>(entropy_between_neighbors[i]);

			bin_entropy[left] = std::get<3>(entropy_between_neighbors[i]);
			bin_entropy[right] = std::get<3>(entropy_between_neighbors[i]);

			//Update neighbors
			size_t l = entropy_between_neighbors.size();
			size_t r = entropy_between_neighbors.size();
			for(size_t j=i;j < entropy_between_neighbors.size();++j) {
				if(left == std::get<1>(entropy_between_neighbors[j])) {
					l = j;
					std::vector<size_t>& nbins = std::get<2>(entropy_between_neighbors[j]);
					nbins.insert(nbins.end(), std::get<2>(entropy_between_neighbors[i]).begin(), std::get<2>(entropy_between_neighbors[i]).end());
					std::sort(nbins.begin(), nbins.end());
					auto last = std::unique(nbins.begin(), nbins.end());
					nbins.erase(last, nbins.end());
					value_type z_n = boost::numeric_cast<value_type>(nbins.size()) / boost::numeric_cast<value_type>(X.size());
					std::get<3>(entropy_between_neighbors[j]) = z_n * conditional_cummulative_entropy_for_sorted_vectors(X, nbins);
					std::get<3>(entropy_between_neighbors[j]) -= bin_entropy[std::get<0>(entropy_between_neighbors[j])];
					std::get<3>(entropy_between_neighbors[j]) -= bin_entropy[std::get<1>(entropy_between_neighbors[j])];
				} else if(right ==  std::get<0>(entropy_between_neighbors[j])) {
					r = j;
					std::vector<size_t>& nbins = std::get<2>(entropy_between_neighbors[j]);
					nbins.insert(nbins.end(), std::get<2>(entropy_between_neighbors[i]).begin(), std::get<2>(entropy_between_neighbors[i]).end());
					std::sort(nbins.begin(), nbins.end());
					auto last = std::unique(nbins.begin(), nbins.end());
					nbins.erase(last, nbins.end());
					value_type z_n = boost::numeric_cast<value_type>(nbins.size()) / boost::numeric_cast<value_type>(X.size());
					std::get<3>(entropy_between_neighbors[j]) = z_n * conditional_cummulative_entropy_for_sorted_vectors(X, nbins);
					std::get<3>(entropy_between_neighbors[j]) -= bin_entropy[std::get<0>(entropy_between_neighbors[j])];
					std::get<3>(entropy_between_neighbors[j]) -= bin_entropy[std::get<1>(entropy_between_neighbors[j])];
				}
			}

			if (l < entropy_between_neighbors.size()) { 
				reorder_vector_inplace(entropy_between_neighbors, i+1, l, r);
			}
			if (r < entropy_between_neighbors.size()) {
				reorder_vector_inplace(entropy_between_neighbors, i+1, r, entropy_between_neighbors.size());
			}

			entropies.emplace_back(entropy);
		}

		return entropies;
	}

		/**
	* This method calculates the conditional cummulative entropy h(X|Y).
	*
	* @param X Vector of scores.
	* @param Y Vector of vector of scores.
	*
	* @return conditional cummulative entropy h(X|Y)
	*/
	template <typename value_type>
	std::vector<value_type> conditional_cummulative_entropy_kruskall_estimator_for_sorted_vectors(std::vector<value_type>& X, std::vector<std::vector<value_type> >& Y)
	{	
		// Compute paiwise distances
		std::vector<std::tuple<size_t, size_t, value_type, value_type> > distances = compute_pairwise_distances(X, Y, SquaredDistance());
		std::vector<std::vector<size_t> > bins(Y.size());
		std::vector<std::vector<size_t>* > bin_ptr(Y.size());
		for(size_t i=0; i<bins.size(); ++i) {
			//Maximal size?
			bins[i].reserve(Y.size()-i);
			//Indices of Xs
			bins[i].emplace_back(i);
			bin_ptr[i] = &bins[i];
		}

		//Kruskall
		std::vector<value_type> entropies;
		entropies.reserve(X.size());
		//Do not forget the first step
		entropies.emplace_back(value_type(0));
		for(size_t i=0; i<distances.size(); ++i){
			size_t left = std::get<0>(distances[i]);
			size_t right = std::get<1>(distances[i]);
			if(bin_ptr[left] == bin_ptr[right] || bin_ptr[left]->size() == 0 || bin_ptr[right]->size() == 0) {
				continue;
			}
			bin_ptr[left]->insert(bin_ptr[left]->end(), bin_ptr[right]->begin(), bin_ptr[right]->end());
			bin_ptr[right]->clear();
			//This might actually be a problem
			std::sort(bin_ptr[left]->begin(), bin_ptr[left]->end());
			
			bin_ptr[right]=bin_ptr[left];
			entropies.emplace_back(conditional_cummulative_entropy_for_sorted_vectors(X, bins));
		}

		return entropies;
	}

	/**
	* This method calculates the conditional cummulative entropy h(X|Y).
	*
	* @param X Vector of scores.
	* @param Y Vector of vector of scores.
	*
	* @return conditional cummulative entropy h(X|Y)
	*/
	template <typename value_type>
	std::vector<value_type> conditional_cummulative_entropy_estimator(const std::vector<value_type>& X, const std::vector<std::vector<value_type> >& Y)
	{
		//Sort vector X
		std::vector<value_type> Xs(std::distance(X.begin(),X.end()));
		std::copy(X.begin(), X.end(), Xs.begin());
		std::sort(Xs.begin(), Xs.end());
		
		//Sort indices according to X
		std::vector<size_t> indices(X.size());	
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b){return X[a] < X[b];});

		// Combine all Ys into one vector and sort this vector with respect to X.
		std::vector<std::vector<value_type> > Y_transformed;
		Y_transformed.reserve(Y[0].size());
		for(size_t i : indices){
			std::vector<value_type> Yi(Y.size());
			for(size_t j=0; j<Y.size(); ++j) {
				Yi[j]=Y[j][i]; 
			}
			Y_transformed.emplace_back(std::move(Yi));
		}
		
		return conditional_cummulative_entropy_estimator_for_sorted_vectors(Xs, Y_transformed);
	}

		/**
	* This method calculates the conditional cummulative entropy h(X|Y).
	*
	* @param X Vector of scores.
	* @param Y Vector of vector of scores.
	*
	* @return conditional cummulative entropy h(X|Y)
	*/
	template <typename value_type>
	std::vector<value_type> conditional_cummulative_entropy_kruskall_estimator(const std::vector<value_type>& X, const std::vector<std::vector<value_type> >& Y)
	{
		//Sort vector X
		std::vector<value_type> Xs(std::distance(X.begin(),X.end()));
		std::copy(X.begin(), X.end(), Xs.begin());
		std::sort(Xs.begin(), Xs.end());
		
		//Sort indices according to X
		std::vector<size_t> indices(X.size());	
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b){return X[a] < X[b];});

		// Combine all Ys into one vector and sort this vector with respect to X.
		std::vector<std::vector<value_type> > Y_transformed;
		Y_transformed.reserve(Y[0].size());
		for(size_t i : indices){
			std::vector<value_type> Yi(Y.size());
			for(size_t j=0; j<Y.size(); ++j) {
				Yi[j]=Y[j][i]; 
			}
			Y_transformed.emplace_back(std::move(Yi));
		}
		
		return conditional_cummulative_entropy_kruskall_estimator_for_sorted_vectors(Xs, Y_transformed);
	}
}
}

#endif // GT2_CORE_ENTROPY_H
