%module gtcore
%{
#include <genetrail2/core/Exception.h>
%}

namespace GeneTrail
{
	class IOError : public std::exception
	{
		public:
			IOError(const std::string& error);
			virtual const char* what() const;

		private:
			std::string msg_;
	};

	class NotImplemented : public std::exception
	{
		public:
			NotImplemented(const char* file, int line, const std::string& method);
			const char* what() const;

		private:
			std::string msg_;
	};

	class InvalidIndex : public std::exception
	{
		public:
			InvalidIndex(unsigned int i, unsigned int max);
			const char* what() const;

		private:
			std::string msg_;
	};
}

%typemap(throws, throws="java.lang.IndexOutOfBoundsException") GeneTrail::InvalidIndex {
#ifdef SWIGJAVA
	jclass excep = jenv->FindClass("java/lang/IndexOutOfBoundsException");
	if (excep)
		jenv->ThrowNew(excep, $1.what());

	return $null;
#elif SWIGPYTHON
    PyErr_SetString(PyExc_RuntimeError, $1.what());
    SWIG_fail;
#endif
}
