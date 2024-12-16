#include "graph.hpp"
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

Graph::Graph(std::vector<Tree *> trees, Taxa &subset, std::string weighting) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    // graph = new weight_t**[2];
    // graph[0] = Matrix::new_mat(size);
    // graph[1] = Matrix::new_mat(size);
    parlay::sequence<weight_t***> subgraphs(trees.size());
    if (verbose > "1") count[1] = count[2] = count[3] = 0;
    parlay::parallel_for(0, trees.size(), [&](size_t k) {
        Taxa thread_local_subset(subset);
        // for (int j = 0; j < k; ++j) {
        //     std::unordered_map<index_t, index_t> valid = trees[j]->get_indices();
        //     thread_local_subset.weight_update(valid);
        // }
        std::unordered_map<index_t, index_t> valid = trees[k]->get_indices();
        thread_local_subset.weight_update(valid);
        subgraphs[k] = new weight_t**[2];
        subgraphs[k][0] = Matrix::new_mat(size);
        subgraphs[k][1] = Matrix::new_mat(size);
        weight_t*** temp;
        if (weighting == "f")
            temp = trees[k]->build_graph(thread_local_subset);
        else
            temp = trees[k]->build_wgraph(thread_local_subset);
        for (index_t i = 0; i < size; i ++) {
            for (index_t j = 0; j < size; j ++) {
                if (temp[0][i][j] > 0 || temp[1][i][j] > 0) {
                    index_t i_ = index2index[subset.root_at(i)];
                    index_t j_ = index2index[subset.root_at(j)];
                    subgraphs[k][0][i_][j_] += temp[0][i][j];
                    subgraphs[k][1][i_][j_] += temp[1][i][j];
                }
            }
        }
    });
    auto identity = new weight_t**[2];
    identity[0] = Matrix::new_mat(size);
    identity[1] = Matrix::new_mat(size);
    // REPLACE THIS WITH PARALLEL SUM OR SOMETHING
    graph = parlay::reduce(subgraphs,
    parlay::binary_op([&](weight_t*** g1, weight_t*** g2)->weight_t*** {
        auto graphy = new weight_t**[2];
        graphy[0] = Matrix::new_mat(size);
        graphy[1] = Matrix::new_mat(size);
        for (index_t i = 0; i < size; i ++) {
            for (index_t j = 0; j < size; j ++) {
                if (g1[0][i][j] > 0 || g1[1][i][j] > 0) {
                    graphy[0][i][j] += g1[0][i][j];
                    graphy[1][i][j] += g1[1][i][j];
                }
                if (g2[0][i][j] > 0 || g2[1][i][j] > 0) {
                    graphy[0][i][j] += g2[0][i][j];
                    graphy[1][i][j] += g2[1][i][j];
                }
            }
        }
        return graphy;
    }, identity) 
    );
    // UPDATE THIS WITH PARALLEL MEMORY ALLOCATOR
    for (int i = 0; i < trees.size(); ++i) {
        Matrix::delete_mat(subgraphs[i][0], size);
        Matrix::delete_mat(subgraphs[i][1], size);
        delete [] subgraphs[i];
    }
    // for (Tree *tree : trees) {
    //     std::unordered_map<index_t, index_t> valid = tree->get_indices();
    //     subset.weight_update(valid);
    //     weight_t ***subgraph;
    //     if (weighting == "f")
    //         subgraph = tree->build_graph(subset);
    //     else
    //         subgraph = tree->build_wgraph(subset);
    //     for (index_t i = 0; i < size; i ++) {
    //         for (index_t j = 0; j < size; j ++) {
    //             if (subgraph[0][i][j] > 0 || subgraph[1][i][j] > 0) {
    //                 index_t i_ = index2index[subset.root_at(i)];
    //                 index_t j_ = index2index[subset.root_at(j)];
    //                 graph[0][i_][j_] += subgraph[0][i][j];
    //                 graph[1][i_][j_] += subgraph[1][i][j];
    //             }
    //         }
    //     }
    //     Matrix::delete_mat(subgraph[0], size);
    //     Matrix::delete_mat(subgraph[1], size);
    //     delete [] subgraph;
    // }
    /*
    weight_t **temp_graph = Matrix::new_mat(size);
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            temp_graph[i][j] = graph[0][i][j] + graph[1][i][j];
        }
    }
    std::cout << Matrix::display_mat(temp_graph, size) << std::endl;
    Matrix::delete_mat(temp_graph, size);
    */
    if (verbose > "1") subproblem_csv << ',' << count[1] << ',' << count[2] << ',' << count[3];
}

void Graph::write_good_edges(Dict *dict) {
    index_t i_, j_;
    weight_t normval, nelem;
    
    // Find average value in good edges   
    /*nelem = size * (size - 1) / 2; // 0 + 1 + 2 + ... + n-1 = n(n-1)/2
    normval = 0.0;
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            normval += (graph[0][i][j] / nelem);
        }
    }*/

    // Find minimum value in good edges
    normval = std::numeric_limits<double>::max();
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            if (graph[0][i][j] > 1 && graph[0][i][j] < normval)
                normval = graph[0][i][j];
        }
    }

    std:: cout << "Normalizing good edges by " << normval << std::endl;

    good_edges_txt << size << std::endl;
    for (size_t i = 0; i < size; i++) {
        i_ = index2index[i];
        good_edges_txt << dict->index2label(i);
        for (size_t j = 0; j < size; j++) {
            j_ = index2index[j];
            good_edges_txt << " " << graph[0][i_][j_] / normval;
        }
        good_edges_txt << std::endl;
    }
}

void Graph::write_bad_edges(Dict *dict) {
    index_t i_, j_;
    weight_t normval, nelem;
    
    // Find average value in good edges   
    /*nelem = size * (size - 1) / 2; // 0 + 1 + 2 + ... + n-1 = n(n-1)/2
    normval = 0.0;
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            normval += (graph[1][i][j] / nelem);
        }
    }*/

    // Find minimum value in good edges
    normval = std::numeric_limits<double>::max();
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            if (graph[1][i][j] > 1 && graph[1][i][j] < normval)
                normval = graph[1][i][j];
        }
    }

    std:: cout << "Normalizing bad edges by " << normval << std::endl;

    bad_edges_txt << size << std::endl;
    for (size_t i = 0; i < size; i++) {
        i_ = index2index[i];
        bad_edges_txt << dict->index2label(i);
        for (size_t j = 0; j < size; j++) {
            j_ = index2index[j];
            bad_edges_txt << " " << graph[1][i_][j_] / normval;
        }
        bad_edges_txt << std::endl;
    }
}

Graph::Graph(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(size);
    graph[1] = Matrix::new_mat(size);
    for (auto elem : quartets) {
        index_t *indices = split(elem.first);
        index_t a = index2index[indices[0]], b = index2index[indices[1]], c = index2index[indices[2]], d = index2index[indices[3]];
        weight_t w = elem.second;
        // bad edges
        graph[1][a][b] += w; graph[1][c][d] += w; graph[1][b][a] += w; graph[1][d][c] += w;

        // good edges
        graph[0][a][c] += w; graph[0][a][d] += w; graph[0][b][c] += w; graph[0][b][d] += w;
        graph[0][c][a] += w; graph[0][d][a] += w; graph[0][c][b] += w; graph[0][d][b] += w;
        delete [] indices;
    }
    /*
    weight_t **temp_graph = Matrix::new_mat(size);
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            temp_graph[i][j] = graph[0][i][j] + graph[1][i][j];
        }
    }
    std::cout << Matrix::display_mat(temp_graph, size) << std::endl;
    Matrix::delete_mat(temp_graph, size);
    */
}

Graph::~Graph() {
    Matrix::delete_mat(graph[0], size);
    Matrix::delete_mat(graph[1], size);
    delete[] graph;
}

std::string Graph::to_string() {
    return Matrix::display_mat(graph[0], size) + "\n" + Matrix::display_mat(graph[1], size);
}

weight_t Graph::get_cut(std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit) {
    weight_t positive_weight = -1.0;
    std::vector<index_t> a, b;
    weight_t lower = 0.0, upper = 6.0;
    /*
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            if (graph[1][i][j] == 0) continue;
            weight_t ratio = graph[0][i][j] / graph[1][i][j];
            if (ratio > upper) upper = ratio;
        }
    }
    */
    while (lower + 0.1 < upper) {
        weight_t alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        weight_t weight = sdp_cut(alpha, &a, &b, iter_limit);
        if (weight < 0.001 || a.size() <= 1 || b.size() <= 1) {
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = alpha;
            *A = a;
            *B = b;
        }
    }
    /*
    if (A->size() <= 1 || B->size() <= 1) {
        std::cout << Matrix::display_mat(graph[0], size) << std::endl;
        std::cout << Matrix::display_mat(graph[1], size) << std::endl;
    }
    assert(A->size() > 1 && B->size() > 1);
    */
    return positive_weight;
}

weight_t Graph::sdp_cut(weight_t alpha, std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit) {
    std::vector<Instance::InstanceTuple> input;
    weight_t avg = 0, num = size * (size - 1) / 2;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t temp = (graph[0][i][j] - alpha * graph[1][i][j]) / num;
            if (temp < 0) temp = - temp;
            avg += temp;
        }
    }
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t weight = (graph[0][i][j] - alpha * graph[1][i][j]) / avg;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    QuartetGraphMaxCutCallback mc(iter_limit);
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, &mc);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (index_t i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->push_back(indices[i]);
        else 
            B->push_back(indices[i]);
    }
    return solution.get_weight();
}

QuartetGraphMaxCutCallback::QuartetGraphMaxCutCallback(unsigned long int iter_limit) {
    this->iter_limit = iter_limit;
}

QuartetGraphMaxCutCallback::~QuartetGraphMaxCutCallback() {

}

bool QuartetGraphMaxCutCallback::Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime) {
    return true;
}

bool QuartetGraphMaxCutCallback::Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime, int iter) {
    assert(iter >= 0);
    return iter < this->iter_limit;
}
