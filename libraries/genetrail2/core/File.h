#ifndef GT2_CORE_FILE_H
#define GT2_CORE_FILE_H

#include <string>
#include <memory>

#include <ios>

#include "macros.h"

namespace GeneTrail
{
	// Note: This generates an warning in gcc. This is
	//       a compiler bug and can be ignored.
	enum class GT2_EXPORT FileOpenMode {
		READ,
		WRITE,
		APPEND
	};

	class GT2_EXPORT GenericFile
	{
		public:
		GenericFile(const std::string& path, FileOpenMode mode);

		GenericFile(const GenericFile&) = delete;
		GenericFile& operator=(const GenericFile&) = delete;

		GenericFile(GenericFile&&) = default;
		GenericFile& operator=(GenericFile&&) = default;

		virtual ~GenericFile();

		operator bool() const;
		bool operator!() const;

		bool isReading() const;
		bool isWriting() const;

		FileOpenMode mode() const;

		protected:
		virtual bool isValid_() const;

		// The main stream from which we will read
		std::istream* in_strm_;
		std::ostream* out_strm_;

		private:
		enum class Compressor {
			BZIP2,
			GZIP,
			NONE,
			ZLIB
		};

		FileOpenMode mode_;

		bool destruct_;

		Compressor getCompressor_(const std::string& path) const;

		void openForReading_(const std::string& path);
		void openForWriting_(const std::string& path, std::ios_base::openmode m);
	};

	template <typename ReturnType> class File : public GenericFile
	{
		public:
		using GenericFile::GenericFile;

		virtual ~File() {};

		virtual ReturnType read() = 0;
		virtual bool write(const ReturnType& r) = 0;

		template<typename R>
		friend File<R>& operator<<(File<R>&, const R&);
		template<typename R>
		friend File<R>& operator>>(File<R>&, R&);
	};

	template <typename ReturnType>
	File<ReturnType>& operator<<(File<ReturnType>& f, const ReturnType& r)
	{
		f.write(r);
		return f;
	}

	template <typename ReturnType>
	File<ReturnType>& operator>>(File<ReturnType>& f, ReturnType& r)
	{
		r = f.read();
		return f;
	}
}

#endif // GT2_CORE_FILE_H

