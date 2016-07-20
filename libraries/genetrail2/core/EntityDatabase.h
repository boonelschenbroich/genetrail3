/*
 * GeneTrail2 - An efficient library for interpreting genetic data
 * Copyright (C) 2015 Daniel Stöckel dstoeckel@bioinf.uni-sb.de>
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

#ifndef GT2_ENTITY_DATABASE_H
#define GT2_ENTITY_DATABASE_H

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <vector>

#include "macros.h"

namespace GeneTrail
{
	/**
	 * This class implements the flyweight pattern for entity names
	 * (e.g. gene names, protein names, ...).
	 *
	 * It also provides facilities for converting a range of input strings to
	 * the respective handles.
	 *
	 * @warning Note that the class is currently not thread-safe.
	 */
	class GT2_EXPORT EntityDatabase
	{
		public:
		EntityDatabase() = default;
		EntityDatabase(const EntityDatabase&) = default;
		EntityDatabase(EntityDatabase&&) = default;

		EntityDatabase& operator=(const EntityDatabase&) = default;
		EntityDatabase& operator=(EntityDatabase&&) = default;

		/**
		 * Removes all entities from the database.
		 */
		void clear();

		/**
		 * Return the name of instance i
		 *
		 * @param i A valid id of an instance.
		 * @return The name of the instance.
		 */
		const std::string& name(size_t i) const { return db_[i]; }

		/**
		 * Return the id of an entity. If the entity is not yet known,
		 * a new id will be generated and returned.
		 *
		 * @param name The name of an entity.
		 * @return The (possibly new) id of the entity.
		 */
		size_t index(const std::string& name);

		/**
		 * Return the id of an entity. This is the constant
		 * version, which throws an exception if the entity is not known.
		 *
		 * @param name The name of an entity.
		 * @return The id of the entity.
		 *
		 * @throw UnkownEntry if the entity is not registered with the database.
		 */
		size_t index(const std::string& name) const;

		/**
		 * Overload for index(const std::string&)
		 */
		size_t operator()(const std::string& name) { return index(name); }

		/**
		 * Overload for index(const std::string&) const
		 */
		size_t operator()(const std::string& name) const { return index(name); }

		/**
		 * Overload for name(size_t i) const
		 */
		const std::string& operator()(size_t i) const { return name(i); }

		/**
		 * Helper function that looks up the entries in the range [a,b) and
		 * inserts the name/index into c.
		 *
		 * @param a The start of the input range.
		 * @param b The end of the input range.
		 * @param c The start of the output range.
		 */
		template <typename InputIterator, typename OutputIterator>
		void transform(InputIterator a, InputIterator b, OutputIterator c)
		{
			std::transform(a, b, c, std::ref(*this));
		}

		/**
		 * Helper function that looks up the entries in the range [a,b) and
		 * inserts the name/index into c.
		 *
		 * @param a The start of the input range.
		 * @param b The end of the input range.
		 * @param c The start of the output range.
		 *
		 * @throw UnkownEntry if [a,b) contains std::string instances that
		 *                    cannot be found in the database.
		 */
		template <typename InputIterator, typename OutputIterator>
		void transform(InputIterator a, InputIterator b, OutputIterator c) const
		{
			std::transform(a, b, c, std::cref(*this));
		}

		/**
		 * Helper function that looks up the entries in the input container and
		 * inserts the name/index into c.
		 *
		 * @param container The container with the input data.
		 * @param out The start of the output range.
		 */
		template <typename Container, typename OutputIterator>
		void transform(const Container& container, OutputIterator out)
		{
			std::transform(container.begin(), container.end(), out, std::ref(*this));
		}

		/**
		 * Helper function that looks up the entries in the input container and
		 * inserts the name/index into c.
		 *
		 * @param container The container with the input data.
		 * @param out The start of the output range.
		 *
		 * @throw UnkownEntry if the container contains std::string instances
		 *                    that cannot be found in the database.
		 */
		template <typename Container, typename OutputIterator>
		void transform(const Container& container, OutputIterator out) const
		{
			std::transform(container.begin(), container.end(), out, std::cref(*this));
		}

		private:
		std::unordered_map<std::string, size_t> name_to_index_;
		std::vector<std::string> db_;
	};
}

#endif // GT2_ENTITY_DATABASE_H
