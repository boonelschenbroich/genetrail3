%module gtcore
%{
#include <genetrail2/core/Matrix.h>
%}

namespace GeneTrail
{
	class Matrix
	{
		public:
			typedef double       value_type;
			typedef unsigned int index_type;

			virtual void setRowNames(const std::vector<std::string>& row_names) = 0;
			virtual void setColNames(const std::vector<std::string>& col_names) = 0;

			virtual const std::vector<std::string>& colNames() const = 0;
			virtual const std::vector<std::string>& rowNames() const = 0;

			virtual void setRowName(const std::string& old_name, const std::string& new_name) = 0;
			virtual void setRowName(index_type i, const std::string& new_name) = 0;

			virtual const std::string& rowName(index_type i) const = 0;

			virtual void setColName(const std::string& old_name, const std::string& new_name) = 0;
			virtual void setColName(index_type j, const std::string& new_name) = 0;
			virtual const std::string& colName(index_type j) const = 0;
			virtual index_type rowIndex(const std::string& row) const = 0;
			virtual index_type colIndex(const std::string& col) const = 0;
			virtual bool hasRow(const std::string& name) const = 0;
			virtual bool hasCol(const std::string& name) const = 0;
			virtual index_type cols() const = 0;
			virtual index_type rows() const = 0;

			value_type set(index_type i, index_type j, value_type v) throw(InvalidIndex);
			value_type get(index_type i, index_type j) const throw(InvalidIndex);

			virtual void shuffleRows(const std::vector<index_type>& perm) = 0;
			virtual void shuffleCols(const std::vector<index_type>& perm) = 0;

			virtual void removeRows(const std::vector<index_type>& indices) = 0;
			virtual void removeCols(const std::vector<index_type>& indices) = 0;

			virtual void transpose() = 0;

	};
}
