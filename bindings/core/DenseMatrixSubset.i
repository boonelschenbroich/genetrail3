%module gtcore
%{
#include <genetrail2/core/DenseMatrixSubset.h>

typedef GeneTrail::DenseMatrix::index_type index_type;
%}

namespace GeneTrail
{
	class DenseMatrixSubset : public Matrix
	{
		public:
			typedef std::vector<DenseMatrix::index_type> ISubset;
			typedef std::vector<std::string> SSubset;

			static DenseMatrixSubset createRowSubset(DenseMatrix* mat, ISubset rows);
			static DenseMatrixSubset createRowSubset(DenseMatrix* mat, const SSubset& rows);
			static DenseMatrixSubset createColSubset(DenseMatrix* mat, ISubset cols);
			static DenseMatrixSubset createColSubset(DenseMatrix* mat, const SSubset& cols);

			DenseMatrixSubset(DenseMatrix* mat, ISubset  rows, ISubset  cols);

			virtual const std::string& colName(index_type j) const;
			virtual const std::string& rowName(index_type i) const;

			virtual index_type colIndex(const std::string& col) const;
			virtual index_type rowIndex(const std::string& row) const;

			virtual index_type cols() const;
			virtual index_type rows() const;

			virtual bool hasCol(const std::string& name) const;
			virtual bool hasRow(const std::string& name) const;

			virtual void setColName(index_type j, const std::string& new_name);
			virtual void setColName(const std::string& old_name, const std::string& new_name);
			virtual void setColNames(const std::vector< std::string >& col_names);
			virtual void setRowName(index_type i, const std::string& new_name);
			virtual void setRowName(const std::string& old_name, const std::string& new_name);
			virtual void setRowNames(const std::vector< std::string >& row_names);

			virtual void removeCols(const std::vector< index_type >& indices);
			virtual void removeRows(const std::vector< index_type >& indices);
			virtual void shuffleCols(const std::vector< index_type >& perm);
			virtual void shuffleRows(const std::vector< index_type >& perm);
			virtual void transpose();

			virtual const std::vector< std::string >& colNames() const;
			virtual const std::vector< std::string >& rowNames() const;
	};
}
