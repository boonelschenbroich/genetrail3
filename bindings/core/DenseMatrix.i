%module gtcore
%{
#include "../../libraries/core/src/DenseMatrix.h"

typedef GeneTrail::DenseMatrix::index_type index_type;
%}

typedef unsigned int index_type;

%feature("notabstract") DenseMatrix;

namespace GeneTrail
{
	class DenseMatrix : public AbstractMatrix
	{
		public:
			DenseMatrix(unsigned int rows, unsigned int cols);
			DenseMatrix(std::vector<std::string> rows, std::vector<std::string> cols);
			DenseMatrix(const DenseMatrix&);

			index_type rows();
			index_type cols();
	};
}
