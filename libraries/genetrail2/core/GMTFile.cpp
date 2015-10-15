#include "GMTFile.h"

#include "Exception.h"
#include "EntityDatabase.h"

#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/find_iterator.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <functional>

using namespace boost;

// A UnaryFunction object as required by the boost::make_transform_iterator
// function.
using Range = iterator_range<std::string::const_iterator>;
struct copy_range_f : public std::unary_function<Range, std::string>
{
	std::string operator()(const Range& r) const
	{
		return trim_copy(copy_range<std::string>(r));
	}
};

namespace GeneTrail
{
	GMTFile::GMTFile(const std::string& path, FileOpenMode mode)
		: CategoryFile(path, mode)
	{
		if(mode == FileOpenMode::READ) {
			advanceLine_();
		}
	}

	Category GMTFile::read()
	{
		if(!isValid_() || !isReading()) {
			throw IOError("File is not open for reading");
		}

		auto split_it =
		    make_split_iterator(next_line_, first_finder("\t", is_equal()));
		decltype(split_it) end_it;

		if(split_it == end_it) {
			//TODO: Better exception
			throw IOError("To few arguments in line");
		}

		auto name = copy_range<std::string>(*split_it++);

		if(split_it == end_it) {
			//TODO: Better exception
			throw IOError("To few arguments in line");
		}

		auto url = copy_range<std::string>(*split_it++);

		auto cit = make_transform_iterator(split_it, copy_range_f());
		auto cend = make_transform_iterator(end_it, copy_range_f());

		// Create the category using our newly created transform iterators
		// we move name and reference, as we do not need them any longer
		Category c(EntityDatabase::global.get(), cit, cend);
		c.setName(std::move(name));
		c.setReference(std::move(url));

		// Prefetch the next line from the file.
		// This is needed in order to determine whether there still are
		// relevant categories
		advanceLine_();

		return c;
	}

	bool GMTFile::write(const Category& cat)
	{
		if(!isValid_() || !isWriting()) {
			throw IOError("File is not open for writing");
		}

		(*out_strm_) << cat.name() << '\t' << cat.reference();

		for(const auto& s : cat) {
			(*out_strm_) << '\t' << s;
		}

		(*out_strm_) << '\n';

		return true;
	}

	void GMTFile::advanceLine_()
	{
		do {
			std::getline(*in_strm_, next_line_);
			trim(next_line_);
		} while(isValid_() && next_line_ == "");
	}
}

