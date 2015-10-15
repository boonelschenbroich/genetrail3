#include <genetrail2/core/GMTFile.h>
#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>

#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <boost/program_options.hpp>
#include <boost/math/distributions/fisher_f.hpp>

#include <iostream>
#include <fstream>

using namespace GeneTrail;
namespace bpo = boost::program_options;
namespace bm = boost::math;

using IVec = std::vector<DenseMatrix::index_type>;

void computeCovariance(const DenseMatrix::DMatrix& sdata,
                       const DenseMatrix::DMatrix& cdata, const IVec& indices,
                       DenseMatrix::DMatrix& cov)
{
	const double factor = 1.0 / (sdata.rows() + cdata.rows() - 2);
	for(size_t i = 0; i < indices.size(); ++i) {
		for(size_t j = i; j < indices.size(); ++j) {
			cov(j, i) = cov(i, j) =
			    (cdata.col(indices[i]).dot(cdata.col(indices[j])) +
			     sdata.col(indices[i]).dot(sdata.col(indices[j]))) *
			    factor;
		}
	}
}

void computePrecision(DenseMatrix::DMatrix& cov)
{
	Eigen::SelfAdjointEigenSolver<DenseMatrix::DMatrix> solver(cov);

	Eigen::VectorXd values = solver.eigenvalues();
	for(int i = 0; i < values.rows(); ++i) {
		values[i] = fabs(values[i]) < 0.01 ? 0.0 : (1.0 / values[i]);
	}

	cov = solver.eigenvectors() * values.asDiagonal() * solver.eigenvectors().transpose();
}

double significane(double t2, int n, int p)
{
	const auto df = n - p - 1;
	t2 *= static_cast<double>(df) / static_cast<double>((n - 2) * p);

	bm::fisher_f F(p, df);
	return bm::cdf(bm::complement(F, t2));
}

double computeEnrichment(const DenseMatrix& cdata, const DenseMatrix& sdata,
                         const Eigen::VectorXd& diff_m, const Category& c)
{
	// Some needed typedefs
	const auto n_a = cdata.rows();
	const auto n_b = sdata.rows();
	const auto n = n_a + n_b;

	IVec indices;
	indices.reserve(c.size());

	for(const auto& s : c.names()) {
		auto i = cdata.colIndex(s);

		if(i != std::numeric_limits<decltype(i)>::max()) {
			indices.push_back(i);
		}
	}

	if(indices.empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const auto p = indices.size();
	const int df = n - p - 1;

	if(df <= 0) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	DenseMatrix::DMatrix cov(c.size(), c.size());
	computeCovariance(cdata.matrix(), sdata.matrix(), indices, cov);
	computePrecision(cov);

	double result = 0.0;

	for(size_t i = 0; i < p; ++i) {
		result += cov(i, i) * diff_m[indices[i]] * diff_m[indices[i]];
		for(size_t j = i + 1; j < p; ++j) {
			result += 2.0 * cov(i,j) * diff_m[indices[i]] * diff_m[indices[j]];
		}
	}

	if(result < 0.0) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	result *= static_cast<double>(n_a * n_b) / static_cast<double>(n);

	return significane(result, n, p);
}

int main(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string categories, control, sample;
	double significance;
	desc.add_options()
		("help,h", "Display this message")
		("signficance,t", bpo::value<double>(&significance)->default_value(0.01), "The critical value for rejecting the H0 hypothesis.")
		("categories,g", bpo::value<std::string>(&categories)->required(), "A file containing the categories to be tested.")
		("control,c", bpo::value<std::string>(&control)->required(), "A matrix containing the control group.")
		("sample,s",  bpo::value<std::string>(&sample)->required(), "A matrix containing the sample group.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);
		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	DenseMatrixReader reader;

	// Read the the control group
	std::ifstream file(control);

	if(!file) {
		std::cerr << "Could not open " << control << " for reading." << std::endl;
		return -1;
	}

	DenseMatrix cdata = reader.read(file);
	file.close();

	// Read the the sample group
	file.open(sample);

	if(!file) {
		std::cerr << "Could not open " << sample << " for reading." << std::endl;
		return -1;
	}
	DenseMatrix sdata = reader.read(file);
	file.close();

	// Center the sample matrix
	// Remark: Do not use auto with Eigen expressions, as they may not be fully evaluated.
	Eigen::RowVectorXd cm = cdata.matrix().colwise().mean();
	cdata.matrix().rowwise() -= cm;

	// Center the sample matrix
	Eigen::RowVectorXd sm = sdata.matrix().colwise().mean();
	sdata.matrix().rowwise() -= sm;

	std::vector<std::pair<std::string, double>> results;

	// Compute the enrichment
	GMTFile input(EntityDatabase::global, categories);
	auto category_db = input.read();
	for(const auto& c : category_db) {
		double enr = computeEnrichment(cdata, sdata, cm - sm, c);

		if(isnan(enr)) {
			std::cerr << "WARNING: Could not compute p-value for " << c.name() << std::endl;
		} else {
			results.emplace_back(c.name(), enr);
		}
	}

	std::sort(results.begin(), results.end(),
	          [](const std::pair<std::string, double>& a,
	             const std::pair<std::string, double>& b) {
		return a.second < b.second;
	});

	for(size_t i = 0; i < results.size(); ++i) {
		double adjusted = results[i].second * static_cast<double>(results.size()) / (i + 1);
		if(adjusted > significance) {
			break;
		}

		std::cout << results[i].first << ": " << adjusted << "\n";
	}

	return 0;
}

