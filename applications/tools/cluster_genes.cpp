#include <DenseMatrix.h>
#include <DenseMatrixReader.h>
#include <SparseMatrix.h>

#include <METISClusterer.h>
#include <NeighborhoodBuilder.h>

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

void dumpGraphviz(std::ostream& out, SparseMatrix& graph, const std::vector<int>& grouping, int num_cluster)
{
	out << "graph {" << std::endl;

	for(unsigned int i = 0; i < graph.rows(); ++i) {
		out << graph.rowName(i) << " [style=filled,fillcolor=\"" << (grouping[i] / (float)num_cluster) << " 1.0 1.0\"]\n";
	}

	for(int k = 0; k < graph.matrix().outerSize(); ++k) {
		for (SparseMatrix::SMatrix::InnerIterator it(graph.matrix(), k); it; ++it) {
			if(k < it.row()) {
				continue;
			}

			out << graph.rowName(k) << " -- " << graph.rowName(it.row()) << "\n";
		}
	}

	out << "}\n";
}

void writeGrouping(std::ostream& out, const DenseMatrix& mat, const std::vector<int>& grouping)
{
	for(unsigned int i = 0; i < grouping.size(); ++i) {
		std::cout << mat.rowName(i) << "\t" << grouping[i] << "\n";
	}
}

void writePartitions(std::ostream& out, const DenseMatrix& mat, const std::vector<int>& grouping)
{
	//TODO: Implement (this depends on subsets of a matrix...)
}

// TODO: Output of the partitions
int main(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string infile, outfile, similarity, fixedpoint, graphviz, partfile;
	unsigned int num_cluster, num_neighbors;

	desc.add_options()
		("help,h", "Display this message")
		("in,i",         bpo::value<std::string>(&infile)->required(), "Input file")
		("out,o",        bpo::value<std::string>(&outfile)->default_value("stdout"), "Output file")
		("partition,p",  bpo::value<std::string>(&partfile), "Write the partitioned data into files. The supplied file name should contain a % wildcard which will be replaced with the partition number. (Not implemented yet)")
		("clusters,c",   bpo::value<unsigned int>(&num_cluster)->required(), "The number of clusters that should be computed")
		("neighbors,n",  bpo::value<unsigned int>(&num_neighbors)->required(), "The number of neighbors in the neighborhood graph")
		("similarity,s", bpo::value<std::string>(&similarity)->default_value("pearson"), "The similarity measure that should be used for building the neighborhood. (Not implemented yet)")
		("fixedpoint,f", bpo::value<std::string>(&fixedpoint)->default_value("linear"), "The encoding that should be used for the conversion to fixed point edge weights. (Not implemented yet)")
		("graphviz,g",   bpo::value<std::string>(&graphviz), "Dump the computed neighborhood graph and partition to the specified file.");

	try {
		bpo::store(bpo::command_line_parser(argc, argv).options(desc).run(), vm);

		if(!vm["help"].empty()) {
			desc.print(std::cout);
			return 0;
		}

		bpo::notify(vm);
	} catch(bpo::error& e) {
		std::cerr << "Error: " << e.what() << "\n";
		desc.print(std::cerr);
		return -1;
	}

	DenseMatrixReader reader;

	std::ifstream input(infile);

	if(!input) {
		std::cerr << "Could not open " << infile << " for reading!" << std::endl;
		return -1;
	}

	std::cout << "Loading matrix ..." << std::endl;

	DenseMatrix mat = reader.read(input);

	input.close();

	if(num_neighbors > mat.rows()) {
		std::cerr << "Too many neighbors requested. (Data has " << mat.rows()
		          << " features, requested are " << num_neighbors
		          << " neighbors.) Aborting!" << std::endl;
		return -2;
	}

	if(num_cluster > mat.rows()) {
		std::cerr << "Too many neighbors requested. (Data has " << mat.rows()
		          << " features, requested are " << num_cluster
		          << " clusters.) Aborting!" << std::endl;
		return -2;
	}

	std::cout << "Building neighborhood graph ..." << std::endl;

	NeighborhoodBuilder nbuilder;
	nbuilder.setNumNeighbors(num_neighbors);

	SparseMatrix graph = nbuilder.build(mat);

	METISClusterer clusterer;

	clusterer.setNumClusters(num_cluster);

	std::cout << "Computing clustering ..." << std::endl;
	clusterer.computeClusters(graph);

	if(!vm["graphviz"].empty()) {
		std::ofstream out(graphviz);
		dumpGraphviz(out, graph, clusterer.grouping(), num_cluster);
	}

	if(outfile == "stderr")
	{
		writeGrouping(std::cerr, mat, clusterer.grouping());
	} else if(outfile == "stdout")
	{
		writeGrouping(std::cout, mat, clusterer.grouping());
	} else
	{
		std::ofstream out(outfile);

		if(!out) {
			std::cerr << "Could not open " << outfile << " for writing.";
			return -1;
		}

		writeGrouping(out, mat, clusterer.grouping());
	}

	if(!vm["partition"].empty()) {
		std::ofstream out(partfile);

		if(!out) {
			std::cerr << "Could not open " << partfile << " for writing.";
			return -1;
		}

		writePartitions(out, mat, clusterer.grouping());
	}

	return 0;
}
