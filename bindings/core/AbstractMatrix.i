%module gtcore
%{
#include <genetrail2/core/AbstractMatrix.h>
%}

namespace GeneTrail
{
	class AbstractMatrix : public Matrix
	{
		public:
			AbstractMatrix(index_type rows, index_type cols);
			AbstractMatrix(std::vector<std::string> rows, std::vector<std::string> cols);
			AbstractMatrix(const AbstractMatrix&);

			void setRowNames(const std::vector<std::string>& row_names);
			void setColNames(const std::vector<std::string>& col_names);
			const std::vector<std::string>& colNames() const;
			const std::vector<std::string>& rowNames() const;
			void setRowName(const std::string& old_name, const std::string& new_name);
			void setRowName(index_type i, const std::string& new_name);
			const std::string& rowName(index_type i) const;
			void setColName(const std::string& old_name, const std::string& new_name);
			void setColName(index_type j, const std::string& new_name);
			const std::string& colName(index_type j) const;

			index_type rowIndex(const std::string& row) const;
			index_type colIndex(const std::string& col) const;

			bool hasRow(const std::string& name) const;
			bool hasCol(const std::string& name) const;

			index_type cols() const;
			index_type rows() const;

			void transpose();
			void shuffleRows(const std::vector<index_type>&) = 0;
	};
}
