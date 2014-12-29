#ifndef GT2_CORE_FILE_H
#define GT2_CORE_FILE_H

#include <string>
#include <memory>

#include <ios>

#include "macros.h"

namespace GeneTrail
{
// This is here to fix the annoying and buggy "type attribute ignored" warning
// present in GCC. Refer to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43407
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"

	/**
	 * An enum that specifies for which purpose a file should be opened.
	 * Currently only reading, writing, and appending is supported.
	 */
	enum class GT2_EXPORT FileOpenMode {
		/// Open the file in read mode.
		READ,
		/// Open the file in write mode.
		WRITE,
		/// Open the file in append mode.
		APPEND
	};

#pragma GCC diagnostic pop

	/**
	 * This is a base class that provides some convenience features like
	 * transparent compression/decompression.
	 *
	 * It should not be used directly, but rather as a base class
	 * from which to create custom parsers/serializers/deserializers.
	 *
	 * The class is not copyable, but assignable.
	 */
	class GT2_EXPORT GenericFile
	{
		public:
		/**
		 * Default constructor. Opens the file at the specified path
		 * using the specified OpenMode.
		 *
		 * @param path The file which should be opened.
		 * @param mode Specifies whether the file should be opened for reading or
		 *             writing.
		 */
		GenericFile(const std::string& path, FileOpenMode mode);

		GenericFile(const GenericFile&) = delete;
		GenericFile& operator=(const GenericFile&) = delete;

		GenericFile(GenericFile&&) = default;
		GenericFile& operator=(GenericFile&&) = default;

		/**
		 * Virtual destructor.
		 */
		virtual ~GenericFile();

		/**
		 * Bool conversion operator
		 *
		 * @returns true when the underlying stream is in a good state.
		 */
		operator bool() const;

		/**
		 * Logical not operator
		 *
		 * @returns true when the underlying stream is in a bad state.
		 */
		bool operator!() const;

		/**
		 * Check the current open mode.
		 *
		 * @returns true iff the file was opened using OpenMode::READ
		 */
		bool isReading() const;

		/**
		 * Check the current open mode.
		 *
		 * @returns true iff the file was opened using
		 *          OpenMode::WRITE || OpenMode::APPEND.
		 */
		bool isWriting() const;

		/**
		 * Get the OpenMode with which the file was opened.
		 *
		 * @returns The OpenMode with which the file was opened.
		 */
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

	/**
	 * An abstract base class for (de-)serializers. It provides a default
	 * input/output operator implementation as well as a basic interface
	 * for reading and writing.
	 */
	template <typename ReturnType> class File : public GenericFile
	{
		public:
		using GenericFile::GenericFile;

		virtual ~File() {};

		/**
		 * Read a ReturnType instance from the stream and return it.
		 *
		 * @returns A valid ReturnType instance.
		 */
		virtual ReturnType read() = 0;

		/**
		 * Writes a ReturnType instance to the stream.
		 *
		 * @returns true iff writing was completed sucessfully.
		 */
		virtual bool write(const ReturnType& r) = 0;

		/**
		 * Global output operator.
		 */
		template<typename R>
		friend File<R>& operator<<(File<R>&, const R&);

		/**
		 * Global input operator.
		 */
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

