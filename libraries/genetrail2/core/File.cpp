#include "File.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/stream.hpp>

#include <iostream>
#include <fstream>

namespace bio = boost::iostreams;

namespace GeneTrail
{
	GenericFile::GenericFile(const std::string& path, FileOpenMode mode)
	    : in_strm_(nullptr), out_strm_(nullptr), mode_(mode), destruct_(false)
	{
		switch(mode) {
			case FileOpenMode::READ:
				openForReading_(path);
				break;
			case FileOpenMode::WRITE:
				openForWriting_(path, std::ios_base::out);
				break;
			case FileOpenMode::APPEND:
				openForWriting_(path, std::ios_base::out | std::ios_base::app);
				break;
		}
	}

	GenericFile::~GenericFile()
	{
		if(destruct_) {
			delete in_strm_;
			delete out_strm_;
		}
	}

	void GenericFile::openForReading_(const std::string& path)
	{
		if(path == "stdin") {
			in_strm_ = &std::cin;
			return;
		}

		destruct_ = true;

		const Compressor compressor = getCompressor_(path);

		if(compressor == Compressor::NONE) {
			in_strm_ = new std::ifstream(path);
			return;
		}

		auto bstream = new bio::filtering_istream();
		in_strm_ = bstream;

		switch(compressor) {
			case Compressor::GZIP:
				bstream->push(bio::gzip_decompressor());
				break;
			case Compressor::BZIP2:
				bstream->push(bio::bzip2_decompressor());
				break;
			case Compressor::ZLIB:
				bstream->push(bio::zlib_decompressor());
				break;
			case Compressor::NONE:
				// This should never be reached
				assert(false);
		}

		bstream->push(bio::file_source(path));
	}

	void GenericFile::openForWriting_(const std::string& path,
	                                  std::ios_base::openmode m)
	{
		if(path == "stdout") {
			out_strm_ = &std::cout;
			return;
		} else if(path == "stderr") {
			out_strm_ = &std::cerr;
			return;
		}

		destruct_ = true;

		Compressor compressor = getCompressor_(path);

		out_strm_ = new std::ofstream(path, m);

		if(compressor == Compressor::NONE) {
			return;
		}

		auto bstream = new bio::filtering_ostream();
		out_strm_ = bstream;

		switch(compressor) {
			case Compressor::GZIP:
				bstream->push(bio::gzip_compressor());
				break;
			case Compressor::BZIP2:
				bstream->push(bio::bzip2_compressor());
				break;
			case Compressor::ZLIB:
				bstream->push(bio::zlib_compressor());
				break;
			case Compressor::NONE:
				// This should never be reached
				assert(false);
		}

		bstream->push(bio::file_sink(path, m));
	}

	GenericFile::Compressor
	GenericFile::getCompressor_(const std::string& path) const
	{
		if(boost::iends_with(path, ".gz")) {
			return Compressor::GZIP;
		} else if(boost::iends_with(path, ".bz2")) {
			return Compressor::BZIP2;
		} else if(boost::iends_with(path, ".zlib")) {
			return Compressor::ZLIB;
		}

		return Compressor::NONE;
	}

	GenericFile::operator bool() const
	{
		return isValid_();
	}

	bool GenericFile::operator!() const
	{
		return !isValid_();
	}

	bool GenericFile::isReading() const
	{
		return mode_ == FileOpenMode::READ;
	}

	bool GenericFile::isWriting() const
	{
		return mode_ == FileOpenMode::WRITE || mode_ == FileOpenMode::APPEND;
	}

	FileOpenMode GenericFile::mode() const
	{
		return mode_;
	}

	bool GenericFile::isValid_() const
	{
		return (in_strm_ && *in_strm_) || (out_strm_ && *out_strm_);
	}
}

