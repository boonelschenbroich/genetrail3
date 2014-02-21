#include "Pathfinder.h"
#include <algorithm>
#include <map>

using namespace GeneTrail;

// If true the matrix in every step is printed
bool debug = true;

void Pathfinder::printMatrix(const std::vector<std::vector<int> >& m) {
    for (std::vector<int> row : m) {
        for (int v : row) {
            std::string s = (v < 0) ? "" : " ";
            std::cout << s << v << ", ";
        }
        std::cout << std::endl;
    }
}

int Pathfinder::computeRunningSum(int bpi, int n, int i, int l) {
    return bpi * n - i*l;
}

bool geneListValid(const std::vector<std::string>& sorted_gene_list) {
    std::set<std::string> sorted;

    for (const auto& s : sorted_gene_list) {
        sorted.insert(s);
    }

    return sorted.size() == sorted_gene_list.size();
}

void Pathfinder::initializeFields(const GraphType& graph, const std::vector<std::string>& sorted_gene_list, int& length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string, std::string>& regulations) {
    assert(geneListValid(sorted_gene_list));
    assert(boost::num_vertices(graph) == sorted_gene_list.size());

    //The number of genes in the list
    const int numberOfGeneIds = sorted_gene_list.size();

    //Map from name to rank in gene list
    for (int i = 0; i < numberOfGeneIds; ++i) {
        name2rank[sorted_gene_list[i]] = i;
    }

    nodes.clear();
    nodes.resize(boost::num_vertices(graph));

    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; vi++) {
        vertex_descriptor vd = *vi;

        auto vertex_id = boost::get(vertex_identifier, graph, vd);

        nodes[name2rank[vertex_id]] = vd;
        vertex_map[vd] = name2rank[vertex_id];
    }

    // Initialize first layer
    std::vector < std::vector<int >> M(std::vector<std::vector<int> > (numberOfGeneIds, std::vector<int> (numberOfGeneIds, 0)));
    M_1 = M;

    // Fill the first layer
    // This is very simple as the gene list is sorted and we have a bijective mapping
    for (int k = 0; k < numberOfGeneIds; ++k) {
        for (int i = k; i < numberOfGeneIds; ++i) {
            M_1[k][i] = 1;
        }
    }

    // Only needed for debugging
    if (debug) {
        printMatrix(M_1);
    }

    // Compute RS for all paths of length 1
    std::vector < std::vector<int >> rs(std::vector<std::vector<int> > (length, std::vector<int> (numberOfGeneIds, 0)));
    running_sums = rs;

    for (int i = 0; i < numberOfGeneIds; ++i) {
        running_sums[0][i] = computeRunningSum(1, numberOfGeneIds, i + 1, 1);
    }

    std::cout << "Layer 1" << std::endl;
}

void Pathfinder::findBestPredecessor(const GraphType& graph, int& best_pred_k, int& best_pred_k_running_sum, int& pred_flag, int l, vertex_descriptor& vd) {
    for (; iei != iei_end; iei++) {
        vertex_descriptor vd_source = source(*iei, graph);
        int source_kv = vertex_map[vd_source];
        int tmp_rs = running_sums[l - 2][source_kv];

        if (pred_flag == 0 || tmp_rs > best_pred_k_running_sum) {
            int kv = vertex_map[vd];

            // Check if k is already on the path
            // We have to avoid cycles
            if (((kv == 0 && M_1[source_kv][kv] == 0) || (M_1[source_kv][kv - 1] == M_1[source_kv][kv])) && M_1[source_kv][kv] != -1) {
                best_pred_k = source_kv;
                best_pred_k_running_sum = tmp_rs;
                ++pred_flag;
            }
        }
    }
}

void Pathfinder::fillNextLayer(int best_pred_k, int k) {
    const int numberOfGeneIds = nodes.size();
    for (int i = 0; i < numberOfGeneIds; ++i) {
        if (k <= i) {
            M_2[k][i] = M_1[best_pred_k][i] + 1;
        } else {
            M_2[k][i] = M_1[best_pred_k][i];
        }
    }
}

void Pathfinder::computeRunningSum(int k, int l, int& max_runnig_sum_k) {
    const int numberOfGeneIds = nodes.size();
    int bpi = 1;

    if (M_2[k][0] == 1) {
        max_runnig_sum_k = computeRunningSum(bpi, numberOfGeneIds, 1, l);
        ++bpi;
    }

    for (int i = 0; i < numberOfGeneIds - 1; ++i) {
        if (M_2[k][i] < M_2[k][i + 1]) {
            int running_sum_k = computeRunningSum(bpi, numberOfGeneIds, i + 2, l);

            if (bpi == 1 || running_sum_k > max_runnig_sum_k) {
                max_runnig_sum_k = running_sum_k;
            }

            ++bpi;
        }
    }
}

void Pathfinder::computeDeregulatedPath(const GraphType& graph, const std::vector<std::string>& sorted_gene_list, int length, std::vector<std::vector<std::string> >& best_paths, std::map<std::string, std::string>& regulations) {

    initializeFields(graph, sorted_gene_list, length, best_paths, regulations);

    // Initialize fields and first layer of the matrix	
    //initializeFields(graph, sorted_gene_list, length, best_paths, regulations);
    const int numberOfGeneIds = nodes.size();
    //Extend the path
    //Fill the layers 2..(length)
    for (int l = 2; l <= length; ++l) {
        int bestk = -1;
        int bestk_running_sum = -1;

        // Contains a mapping for each vertex to the best predecessor
        std::map<int, int> best_preds;

        // Initialize second layer
        std::vector < std::vector<int >> M2(std::vector<std::vector<int> > (numberOfGeneIds, std::vector<int> (numberOfGeneIds, -1)));
        M_2 = M2;

        for (int k = 0; k < numberOfGeneIds; ++k) {
            // Find vertex that corresponds to k
            vertex_descriptor vd = nodes[k];

            // Get all ingoing edges
            tie(iei, iei_end) = boost::in_edges(vd, graph);

            // Continue if there are no ingoing edges
            if (iei == iei_end) {
                continue;
            }

            // For each target
            // Hold the best values for predecessor
            int best_pred_k;
            int best_pred_k_running_sum;
            int pred_flag = 0;

            // Find the best predecessor
            findBestPredecessor(graph, best_pred_k, best_pred_k_running_sum, pred_flag, l, vd);

            // In case there are only cycles possible
            if (pred_flag == 0) {
                //std::cout << "only cycles possible" << std::endl;
                continue;
            }

            // Save for each k the best predecessor 
            best_preds[k] = best_pred_k;

            if (debug) {
                std::cout << sorted_gene_list[best_pred_k] << "-->" << sorted_gene_list[k] << std::endl;
            }

            // Fill the next layer
            fillNextLayer(best_pred_k, k);


            // Compute the running sum for the current path
            int max_runnig_sum_k;
            computeRunningSum(k, l, max_runnig_sum_k);

            // Save the best running sum
            running_sums[l - 1][k] = max_runnig_sum_k;

            if (max_runnig_sum_k > bestk_running_sum) {
                bestk = k;
                bestk_running_sum = max_runnig_sum_k;
            }
        }

        // Print all running sums
        if (debug) {
            printMatrix(running_sums);
        }

        if (bestk == -1) {
            std::cout << "WARNING: No path of length " << l << " found" << std::endl;
            break;
        }

        best_preds_v.push_back(best_preds);

        if (debug) {
            std::cout << std::endl;
            printMatrix(M_2);
        }

        // Estimate the reversed order
        std::vector<std::string> path;
        std::vector<std::string> rev_path;
        rev_path.push_back(sorted_gene_list[bestk]);

        int tmp_k = bestk;
        for (int i = l - 2; i >= 0; --i) {
            best_preds = best_preds_v[i];
            tmp_k = best_preds[tmp_k];
            rev_path.push_back(sorted_gene_list[tmp_k]);
        }

        // Revert the path
        // std::vector<std::string> path;
        if (debug) {
            std::cout << "Best path: ";
        }
        for (int i = rev_path.size() - 1; i >= 0; --i) {
            path.push_back(rev_path[i]);
            if (debug) {
                std::cout << rev_path[i] << "-->";
            }
        }
        if (debug) {
            std::cout << std::endl;
        }

        // Save all edges with regulations
        for (int i = 0; i < (signed)path.size() - 1; ++i) {
            // get edge regulation
            bool exists;
            vertex_descriptor vd1 = nodes[name2rank[path[i]]];
            vertex_descriptor vd2 = nodes[name2rank[path[i + 1]]];

            edge_descriptor ed;
            boost::tie(ed, exists) = boost::edge(vd1, vd2, graph);
            if (exists) {
                regulations[path[i] + path[i + 1]] = boost::get(edge_regulation_type, graph, ed);
            } else {
                std::cout << "ERROR: No edge between " << path[i] << " " << boost::get(vertex_identifier, graph, vd1) << " and " << path[i + 1] << " " << boost::get(vertex_identifier, graph, vd2) << std::endl;
            }
        }

        best_paths.push_back(path);

        // As each layer only depends on the layer before we only need to save two layers at once
        std::swap(M_1, M_2);

        std::cout << "Layer " << l << std::endl;
    }
}
