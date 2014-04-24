#include <genetrail2/core/DenseMatrix.h>
#include <genetrail2/core/DenseMatrixReader.h>
#include <genetrail2/core/DenseMatrixSubset.h>
#include <genetrail2/core/DenseMatrixWriter.h>
#include <genetrail2/core/SparseMatrix.h>
#include <genetrail2/core/SparseMatrixReader.h>
#include <genetrail2/core/SparseMatrixWriter.h>

#include <genetrail2/cluster/METISClusterer.h>
#include <genetrail2/cluster/NeighborhoodBuilder.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using namespace GeneTrail;
namespace bpo = boost::program_options;

void dumpGraphviz(std::ostream& out, SparseMatrix& graph, const std::vector<int>& grouping, int num_cluster)
{
	out << "graph {" << std::endl;

	for(int i = 0; i < num_cluster; ++i)
	{
		out << "subgraph cluster_" << boost::lexical_cast<std::string>(i) << " {\n";
		for(int k = 0; k < graph.matrix().outerSize(); ++k) {
			if(grouping[k] != i) {
				continue;
			}

			out << "\"" << graph.rowName(k) << "\" [style=filled,fillcolor=\"" << (grouping[k] / (float)num_cluster) << " 1.0 1.0\"]\n";

			for (SparseMatrix::SMatrix::InnerIterator it(graph.matrix(), k); it; ++it) {
				if(k < it.row()) {
					continue;
				}

				if(grouping[it.row()] != i) {
					continue;
				}

				out << "\"" << graph.rowName(k) << "\" -- " << "\"" <<graph.rowName(it.row()) << "\"\n";
			}
		}
		out << "}\n";
	}

	for(int k = 0; k < graph.matrix().outerSize(); ++k) {
		for (SparseMatrix::SMatrix::InnerIterator it(graph.matrix(), k); it; ++it) {
			if(k < it.row()) {
				continue;
			}

			if(grouping[k] == grouping[it.row()]) {
				continue;
			}

			out << "\"" << graph.rowName(k) << "\" -- " << "\"" <<graph.rowName(it.row()) << "\"\n";
		}
	}

	out << "}\n";
}

void writeGrouping(std::ostream& out, const DenseMatrix& mat, const std::vector<int>& grouping)
{
	for(unsigned int i = 0; i < grouping.size(); ++i) {
		out << mat.rowName(i) << "\t" << grouping[i] << "\n";
	}
}

void writePartitions(const std::string& partfile, DenseMatrix& mat, const std::vector<int>& grouping, int num_cluster)
{
	DenseMatrixSubset::ISubset subset;
	subset.reserve(mat.rows() / num_cluster);

	int npad = 0;
	int tmp = num_cluster - 1;

	for(int cmp = 1; cmp < tmp; cmp *= 10, ++npad) {}

	for(int i = 0; i < (int)num_cluster; ++i) {
		subset.resize(0);
		int k = 0;
		for(auto j : grouping) {
			if(j == i) {
				subset.push_back(k);
			}
			++k;
		}

		std::stringstream ss;
		ss << std::setw(npad) << std::setfill('0') << i;

		std::ofstream outfile(boost::replace_first_copy(partfile, "%", ss.str()));

		DenseMatrixSubset sub = DenseMatrixSubset::createRowSubset(&mat, subset);

		DenseMatrixWriter writer;
		writer.writeBinary(outfile, sub);
	}
}

// TODO: Output of the partitions
int main(int argc, char* argv[])
{
	bpo::variables_map vm;
	bpo::options_description desc;

	std::string infile, outfile, similarity, fixedpoint, graphviz, partfile, save_graph, load_graph;
	unsigned int num_cluster, num_neighbors;

	desc.add_options()
		("help,h", "Display this message")
		("in,i",         bpo::value<std::string>(&infile)->required(), "Input file")
		("out,o",        bpo::value<std::string>(&outfile)->default_value("stdout"), "Output file")
		("partition,p",  bpo::value<std::string>(&partfile), "Write the partitioned data into files. The supplied file name should contain a % wildcard which will be replaced with the partition number.")
		("clusters,c",   bpo::value<unsigned int>(&num_cluster)->required(), "The number of clusters that should be computed")
		("neighbors,n",  bpo::value<unsigned int>(&num_neighbors)->required(), "The number of neighbors in the neighborhood graph")
		("similarity,s", bpo::value<std::string>(&similarity)->default_value("pearson"), "The similarity measure that should be used for building the neighborhood. (Not implemented yet)")
		("fixedpoint,f", bpo::value<std::string>(&fixedpoint)->default_value("linear"), "The encoding that should be used for the conversion to fixed point edge weights. (Not implemented yet)")
		("graphviz,g",   bpo::value<std::string>(&graphviz), "Dump the computed neighborhood graph and partition to the specified file.")
		("print-scores,x", "Print the achieved cluster scores.")
		("load-graph,l", bpo::value<std::string>(&load_graph), "Load the neighborhood graph from file.")
		("save-graph,d", bpo::value<std::string>(&save_graph), "Save the neighborhood graph to a file.");

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

	if(num_cluster > mat.rows()) {
		std::cerr << "Too many neighbors requested. (Data has " << mat.rows()
		<< " features, requested are " << num_cluster
		<< " clusters.) Aborting!" << std::endl;
		return -2;
	}

	SparseMatrix graph(0,0);

	if(vm["load-graph"].empty()) {
		if(num_neighbors > mat.rows()) {
			std::cerr << "Too many neighbors requested. (Data has " << mat.rows()
			<< " features, requested are " << num_neighbors
			<< " neighbors.) Aborting!" << std::endl;
			return -2;
		}

		std::cout << "Building neighborhood graph ..." << std::endl;

		NeighborhoodBuilder nbuilder;
		nbuilder.setNumNeighbors(num_neighbors);

		graph = nbuilder.build(std::move(mat));
	} else {
		std::cout << "Loading neighborhood graph ..." << std::endl;

		std::ifstream input(load_graph);
		SparseMatrixReader reader;

		graph = reader.read(input);
	}

	if(!vm["save-graph"].empty()) {
		SparseMatrixWriter writer;
		std::ofstream out(save_graph);

		writer.writeBinary(out, graph);
	}

	graph.matrix().makeCompressed();

	std::cout << "Computing clustering ..." << std::endl;
	METISClusterer clusterer;
	clusterer.setAlgorithm(METISClusterer::Recursive);
	clusterer.setNumClusters(num_cluster);
	clusterer.computeClusters(graph);

	if(!vm["print-scores"].empty()) {
		std::vector<double> intra_scores(num_cluster);
		std::vector<double> inter_scores(num_cluster);
		std::vector<int> num_nodes(num_cluster);

		for(int k = 0; k < graph.matrix().outerSize(); ++k) {
			++num_nodes[clusterer.grouping()[k]];
			for(SparseMatrix::SMatrix::InnerIterator it(graph.matrix(), k); it; ++it) {
				if(clusterer.grouping()[it.row()] != clusterer.grouping()[it.col()]) {
					intra_scores[clusterer.grouping()[it.row()]] += it.value();
				} else {
					inter_scores[clusterer.grouping()[it.row()]] += it.value();
					inter_scores[clusterer.grouping()[it.col()]] += it.value();
				}
			}
		}

		for(int i = 0; i < num_cluster; ++i) {
			std::cout << "Cluster " << i << "\tIntra: " << intra_scores[i] << "\tInter: " << inter_scores[i] << "\tSize: " << num_nodes[i] << std::endl;
		}
	}

	if(!vm["graphviz"].empty()) {
		std::ofstream out(graphviz);
		dumpGraphviz(out, graph, clusterer.grouping(), num_cluster);
	}

	if(outfile == "stderr") {
		writeGrouping(std::cerr, mat, clusterer.grouping());
	} else if(outfile == "stdout") {
		writeGrouping(std::cout, mat, clusterer.grouping());
	} else {
		std::ofstream out(outfile);

		if(!out) {
			std::cerr << "Could not open " << outfile << " for writing.";
			return -1;
		}

		writeGrouping(out, mat, clusterer.grouping());
	}

	if(!vm["partition"].empty()) {
		writePartitions(partfile, mat, clusterer.grouping(), num_cluster);
	}

	return 0;
}
