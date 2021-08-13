
/******************************************************************************
 * vcc.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <utility> // for std::pair
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <climits>

#include <boost/functional/hash.hpp>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "mis/initial_mis/greedy_mis.h"

#include "mis/kernel/branch_and_reduce_algorithm.h"

#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_vcc/redu_vcc.h"
#include "redu_vcc/reducer.h"
#include "branch_and_reduce/b_and_r.h"

#include "sigmod_mis/Graph.h"
//
// #include "mis/ils/ils.h"
#include "mis/mis_config.h"

// #include "redu_vcc/graph_vcc.h"

namespace std
{
    template<> struct hash<std::pair<NodeID,NodeID>>
    {
        std::size_t operator()(std::pair<NodeID,NodeID> const & node_pair) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, std::get<0>(node_pair));
            boost::hash_combine(seed, std::get<1>(node_pair));
            return seed;
        }
    };
};

void transform_graph(graph_access const &G, 
                     graph_access       &transformed_graph,
                     std::vector<std::pair<NodeID,NodeID>> &vertex_to_edge) {

        // Step 1: acyclic orientation
        std::vector<NodeID> order;
        order.reserve(G.number_of_nodes());
        for (NodeID node = 0; node < G.number_of_nodes(); node++) {
            order.push_back(node);
        }

        // Step 2: line graph
        std::unordered_set<std::pair<NodeID,NodeID>> neighbor_hash(G.number_of_edges());
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                neighbor_hash.insert(std::pair<NodeID,NodeID>(u,v));
                neighbor_hash.insert(std::pair<NodeID,NodeID>(v,u));
            } endfor
        } endfor

        vertex_to_edge.clear();
        vertex_to_edge.reserve(G.number_of_edges() / 2);
        std::unordered_map<std::pair<NodeID, NodeID>, NodeID> edge_to_vertex(G.number_of_edges() / 2);
        std::vector<std::vector<NodeID>> edge_to_edges(G.number_of_edges() / 2);

        std::vector<NodeID> const empty;

        std::cout << "Making transformed graph..." << std::endl;
        std::size_t full_num_edges = 0;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                std::pair<NodeID,NodeID> node_pair(u, v);
                if (order[v] < order[u]) {
                    node_pair = std::pair<NodeID,NodeID>(v, u);
                }
                if (edge_to_vertex.find(node_pair) == edge_to_vertex.end()) {
                    vertex_to_edge.push_back(node_pair);
                    edge_to_vertex[node_pair] = vertex_to_edge.size() - 1;
                    //edge_to_edges.push_back(empty);
                }
            } endfor

            forall_out_edges(G, e1, v) {
                NodeID u = G.getEdgeTarget(e1);
                std::pair<NodeID,NodeID> first_pair(v, u);
                if (order[u] < order[v]) {
                    first_pair = std::pair<NodeID,NodeID>(u, v);
                }

                forall_out_edges(G, e2, v) {
                    NodeID w = G.getEdgeTarget(e2);
                    if (w < u) continue;
                    if (w == u) continue;
                    std::pair<NodeID,NodeID> second_pair(v, w);
                    if (order[w] < order[v]) {
                        second_pair = std::pair<NodeID,NodeID>(w, v);
                    }

                    bool bad = std::get<0>(first_pair) == std::get<0>(second_pair);
                    std::pair<NodeID,NodeID> const bad_edge = 
                        std::pair<NodeID,NodeID>(std::get<1>(first_pair), std::get<1>(second_pair));
                    bad = bad and (neighbor_hash.find(bad_edge) != neighbor_hash.end());
                    if (bad) continue;

                    // Step 3: filter triples, TODO/DS: appears to work
                    edge_to_edges[edge_to_vertex[first_pair]].push_back(edge_to_vertex[second_pair]);
                    edge_to_edges[edge_to_vertex[second_pair]].push_back(edge_to_vertex[first_pair]);
                    full_num_edges++;

                } endfor
            } endfor
        } endfor

        // Step 3: trim triples (integrated above)

        std::cout << "Sorting neighborhoods..." << std::endl;
        std::size_t num_edges = 0;
        for (std::vector<NodeID> &neighbors : edge_to_edges) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        std::cout << "Transformed Graph has " << edge_to_edges.size() << " vertices and " << num_edges << " edges" << std::endl;

        std::cout << "Constructing graph_access data structure..." << std::endl;

        transformed_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = transformed_graph.new_node();
            transformed_graph.setPartitionIndex(shadow_node, 0);
            transformed_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = transformed_graph.new_edge(shadow_node, neighbor);
                transformed_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        transformed_graph.finish_construction();

        //graph_io::writeGraph(transformed_graph, "transformed.graph");

        //graph_access new_transformed_graph;
        //graph_io::readGraphWeighted(new_transformed_graph, "transformed-sorted.graph");
}


void transform_is_to_clique_cover(
        graph_access &transformed_graph, 
        std::vector<std::pair<NodeID, NodeID>> const &vertex_to_edge,
        std::size_t const num_vertices_original_graph,
        NodeID *solution_is,
        std::size_t const num_cliques,
        std::vector<std::vector<NodeID>> &clique_cover) {

    std::cout << "Reconstructing clique cover..." << std::endl;
    clique_cover.clear();
    clique_cover.reserve(num_cliques);
    unordered_map<NodeID, NodeID> vertex_to_clique_id(num_cliques);
    std::vector<bool> covered(num_vertices_original_graph, false);

    forall_nodes(transformed_graph, v) {
        if (solution_is[v] == 1) {
            std::pair<NodeID,NodeID> const &edge = vertex_to_edge[v];
            //std::cout << "Edge " << std::get<0>(edge) << "," << std::get<1>(edge)
            //         << " is in the independent set" << std::endl;
            if (vertex_to_clique_id.find(std::get<0>(edge)) 
                    == vertex_to_clique_id.end()) {
                vertex_to_clique_id[std::get<0>(edge)] = clique_cover.size();
                clique_cover.push_back(std::vector<NodeID>{std::get<0>(edge)});
                covered[std::get<0>(edge)] = true;
            }
            clique_cover[vertex_to_clique_id[std::get<0>(edge)]].push_back(std::get<1>(edge));
            covered[std::get<1>(edge)] = true;
        }
    } endfor

    for (NodeID v = 0; v < num_vertices_original_graph; v++) {
        if (!covered[v]) {
            //std::cout << "DS: Vertex " << v << " is uncovered...covering with singleton" << std::endl;
            clique_cover.emplace_back(std::vector<NodeID>{v});
        }
    }
    std::cout << "clique_cover has size: " << clique_cover.size() << std::endl;
}

void run_ils(ils &ils_instance, 
             PartitionConfig const &partition_config,
             graph_access &graph, 
             timer &the_timer,
             std::size_t const is_offset) {
        std::cout << "Performing ILS..." << std::endl;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils_instance.perform_ils(ils_config, graph, ils_config.ils_iterations, 
                                 the_timer.elapsed(), is_offset, partition_config.mis);
}

void run_peeling(graph_access &graph, 
                 graph_access &peeled_graph,
                 std::vector<NodeID> &new_to_old_id,
                 std::size_t &cover_offset) {
    Graph mis_G;
    mis_G.read_graph(graph);
    std::vector<std::vector<NodeID>> kernel;
    cover_offset = 0;
    std::cout << "Running peeling..." << endl;
    std::vector<bool> in_initial_is;
    unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id, in_initial_is, cover_offset);
    std::cout << "    Reduced " << graph.number_of_nodes() << " -> " << kernel.size() << " nodes" << std::endl;
    std::cout << "    Offset =  " << cover_offset << std::endl;

    std::size_t num_edges = 0;
    for (std::vector<NodeID> & neighbors : kernel) {
        std::sort(neighbors.begin(), neighbors.end());
        num_edges += neighbors.size();
    }
    num_edges /= 2;

    peeled_graph.set_partition_count(2);
    peeled_graph.start_construction(kernel.size(), 2 * num_edges);
    for (NodeID v = 0; v < kernel.size(); v++) {
        NodeID shadow_node = peeled_graph.new_node();
        peeled_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
        peeled_graph.setNodeWeight(shadow_node, 1);
        for (NodeID const neighbor : kernel[v]) {
            EdgeID shadow_edge = peeled_graph.new_edge(shadow_node, neighbor);
            peeled_graph.setEdgeWeight(shadow_edge, 1);
        }
    }

    peeled_graph.finish_construction();
}


int main(int argn, char **argv) {

    PartitionConfig partition_config;
    std::string graph_filename;

    bool is_graph_weighted = false;
    bool suppress_output   = false;
    bool recursive         = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if(ret_code) {
        return 0;
    }

    std::streambuf* backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if(suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }

    partition_config.LogDump(stdout);
    graph_access G;

    timer t;
    graph_io::readGraphWeighted(G, graph_filename);

    timer s;

    // redu_vcc gVCC(G);
    // std::vector<unsigned int> iso_degree;
    // iso_degree.assign(G.number_of_nodes(), 0);
    // std::vector<unsigned int> dom_degree;
    // dom_degree.assign(G.number_of_nodes(), 0);
    // reducer R;
    // R.exhaustive_reductions(gVCC, iso_degree, dom_degree);
    // gVCC.analyzeGraph(graph_filename, G, s);
    // Graph gra;
    // // std::cout<< "comp"<<std::endl;
    // gra.read_graph(gVCC);
    // // std::cout<< "comp"<<std::endl;
    // unsigned int mis = gra.degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction();
    // std::cout<< mis << std::endl;
    // // gVCC.printAdjList();
    // // std::vector<redu_vcc> comp = gVCC.decompose();
    // // for (redu_vcc &child : comp) {
    // //   child.printAdjList();
    // //   std::cout<<std::endl;
    // // }
    // //
    // return;

    if (partition_config.run_type == "Chalupa") {
        timer total_timer;
        redu_vcc reduVCC(G);
        reduVCC.build_cover();
        double time_to_solution = 0.0;
        reduVCC.solveKernel(partition_config, total_timer, time_to_solution, 0 /* clique cover offset */);
        reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "input_graph_vertices=" << G.number_of_nodes() << std::endl;
        std::cout << "input_graph_edges=" << G.number_of_edges() / 2 << std::endl;
        std::cout << "total_time_to_best=" << time_to_solution << std::endl;
        std::cout << "clique_cover_size=" << reduVCC.clique_cover.size() << std::endl;
        std::cout << "verified_cover=" << (reduVCC.validateCover(G) ? "passed" : "failed") << std::endl;
        std::cout << "optimal=" << (reduVCC.clique_cover.size() == partition_config.mis ? "yes" : "unknown") << std::endl;
        return 0;
    }
    else if (partition_config.run_type == "Redu") {
        redu_vcc reduVCC(G);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(G.number_of_nodes(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(G.number_of_nodes(), 0);
        reducer R;
        R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);
        reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        return 0;
    } else if (partition_config.run_type == "ReduVCC") {
        redu_vcc reduVCC(G);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(G.number_of_nodes(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(G.number_of_nodes(), 0);
        reducer R;
        R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);

        reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);
        reduVCC.build_cover();
        double time_to_solution = 0.0;
        reduVCC.solveKernel(partition_config, s, time_to_solution, R.get_cover_size_offset());
        timer unwind_timer;
        R.unwindReductions(reduVCC, time_to_solution);
        double time_to_unwind = unwind_timer.elapsed();
        reduVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "input_graph_vertices=" << G.number_of_nodes() << std::endl;
        std::cout << "input_graph_edges=" << G.number_of_edges() / 2 << std::endl;
        std::cout << "reduced_graph_vertices=" << reduVCC.kernel_adj_list.size() << std::endl;
        std::cout << "reduced_graph_edges=" << reduVCC.kernel_edges << std::endl;
        std::cout << "total_time_to_best=" << time_to_solution << std::endl;
        std::cout << "time_to_best_without_unwind=" << time_to_solution - time_to_unwind << std::endl;
        std::cout << "clique_cover_size=" << reduVCC.clique_cover.size() << std::endl;
        std::cout << "verified_cover=" << (reduVCC.validateCover(G) ? "passed" : "failed") << std::endl;
        std::cout << "optimal=" << (reduVCC.clique_cover.size() == partition_config.mis ? "yes" : "unknown") << std::endl;

        return 0;
    } else if (partition_config.run_type == "transform") {

        timer transformation_timer;
        // Step 1: acyclic orientation
        std::vector<NodeID> order;
        order.reserve(G.number_of_nodes());
        for (NodeID node = 0; node < G.number_of_nodes(); node++) {
            order.push_back(node);
        }

        // Step 2: line graph
        std::unordered_set<std::pair<NodeID,NodeID>> neighbor_hash;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                neighbor_hash.insert(std::pair<NodeID,NodeID>(u,v));
                neighbor_hash.insert(std::pair<NodeID,NodeID>(v,u));
            } endfor
        } endfor

        std::unordered_map<std::pair<NodeID, NodeID>, NodeID> edge_to_vertex;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        std::vector<std::vector<NodeID>> edge_to_edges;

        std::vector<NodeID> const empty;

        std::cout << "Making filtered line graph..." << std::endl;
        std::size_t full_num_edges = 0;
        forall_nodes(G, v) {
            forall_out_edges(G, e, v) {
                NodeID u = G.getEdgeTarget(e);
                std::pair<NodeID,NodeID> node_pair(u, v);
                if (order[v] < order[u]) {
                    node_pair = std::pair<NodeID,NodeID>(v, u);
                }
                if (edge_to_vertex.find(node_pair) == edge_to_vertex.end()) {
                    vertex_to_edge.push_back(node_pair);
                    edge_to_vertex[node_pair] = vertex_to_edge.size() - 1;
                    edge_to_edges.push_back(empty);
                }
            } endfor

            forall_out_edges(G, e1, v) {
                NodeID u = G.getEdgeTarget(e1);
                std::pair<NodeID,NodeID> first_pair(v, u);
                if (order[u] < order[v]) 
                    first_pair = std::pair<NodeID,NodeID>(u, v);

                forall_out_edges(G, e2, v) {
                    NodeID w = G.getEdgeTarget(e2);
                    if (w <= u) continue;
                    std::pair<NodeID,NodeID> second_pair(v, w);
                    if (order[w] < order[v]) {
                        second_pair = std::pair<NodeID,NodeID>(w, v);
                    }

                    // Step 3: filter triples, TODO/DS: appears to work
                    if (std::get<0>(first_pair) == std::get<0>(second_pair) && 
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(first_pair),std::get<1>(second_pair))) != neighbor_hash.end()) && 
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(second_pair),std::get<1>(first_pair))) != neighbor_hash.end())) continue;

                    
                    edge_to_edges[edge_to_vertex[first_pair]].push_back(edge_to_vertex[second_pair]);
                    edge_to_edges[edge_to_vertex[second_pair]].push_back(edge_to_vertex[first_pair]);
        full_num_edges++;

                } endfor
            } endfor
        } endfor

        // Step 3: trim triples (integrated above)

        double transformation_time = transformation_timer.elapsed();

        // Step 4: Compute MIS
        if (false) // just peeling
        {
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        unsigned int res_mis = mis_G.degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction();

        std::size_t const cover_upper_bound = G.number_of_nodes() - res_mis;
        std::cout << "Original Graph = " << G.number_of_nodes() << ", MIS=" << res_mis << std::endl;
        std::cout << "Transform (Near Linear) Upper Bound = " << cover_upper_bound << std::endl;
        }

        if (false) // just ILS
        {
        std::size_t num_edges = 0;
        for (std::vector<NodeID> &neighbors : edge_to_edges) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        std::cout << "Graph has " << edge_to_edges.size() << " vertices and " << num_edges << " edges" << std::endl;

        std::cout << "Constructing graph_access data structure..." << std::endl;

        graph_access new_graph;
        new_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        //graph_io::writeGraph(new_graph, "transformed.graph");

        //graph_access new_new_graph;
        //graph_io::readGraphWeighted(new_new_graph, "transformed-sorted.graph");

        greedy_mis greedy;
        greedy.initial_partition(partition_config.seed, new_graph);

        std::cout << "Preparing ILS..." << std::endl;
        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::cout << "Performing ILS..." << std::endl;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = G.number_of_nodes() - ils_instance.get_best_solution_size();

        std::cout << "Transform (ILS) Upper Bound = " << ils_cover << std::endl;
        }

        if (false) // just branch-and-reduce
        {
        std::cout << "Building a sorted int version of adjacency list..." << std::endl;
        vector<vector<int>> adjlist(edge_to_edges.size());
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            std::sort(edge_to_edges[v].begin(), edge_to_edges[v].end());
            adjlist[v].reserve(edge_to_edges[v].size());
            for (NodeID const neighbor : edge_to_edges[v])
                adjlist[v].push_back(neighbor);
        }
        branch_and_reduce_algorithm mis_bnr(adjlist, adjlist.size());

        timer bnr_timer;
        bool timeout = mis_bnr.solve(bnr_timer, partition_config.solver_time_limit) == -1;

        std::size_t const exact_cover = G.number_of_nodes() - mis_bnr.get_current_is_size();

        std::cout << "BNR Running time: " << bnr_timer.elapsed() << std::endl;
        std::cout << "Transform (BNR) Upper Bound = " << exact_cover << std::endl;
        std::cout << "    Upper bound is " << (timeout ? "inexact" : "exact") << endl;
        }

        if (true) // peeling + ILS
        {
        timer total_timer;
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        std::vector<std::vector<NodeID>> kernel;
        std::size_t offset = 0;

        std::cout << "Running reducing-peeling..." << endl;
        std::vector<bool> in_initial_is;
        timer peeling_timer;
        std::vector<NodeID> new_to_old_id_unused;
        unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id_unused, in_initial_is, offset);
        double peeling_time = peeling_timer.elapsed();
        std::cout << "Done with reducing-peeling..." << endl;
        std::cout << "    Reduced " << edge_to_edges.size() << " -> " << kernel.size() << " nodes" << std::endl;
        std::cout << "    Offset =  " << offset << std::endl;

        std::size_t num_edges = 0;
        for (std::vector<NodeID> & neighbors : kernel) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        graph_access new_graph;
        new_graph.set_partition_count(2);
        new_graph.start_construction(kernel.size(), 2 * num_edges);
        for (NodeID v = 0; v < kernel.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : kernel[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::size_t solution_offset = G.number_of_nodes() - offset;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations, total_timer.elapsed(), solution_offset, partition_config.mis);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const near_linear_ils_cover = G.number_of_nodes() - (ils_instance.get_best_solution_size() + offset);

        std::cout << "Transform (Near Linear / ILS) Upper Bound = " << near_linear_ils_cover << std::endl;

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "transformation_time=" << transformation_time << std::endl;
        std::cout << "transformed_graph_vertices=" << edge_to_edges.size() << std::endl;
        std::cout << "transformed_graph_edges=" << full_num_edges << std::endl;
        std::cout << "peeling_time=" << peeling_time << std::endl;
        std::cout << "reduced_transformed_graph_vertices=" << new_graph.number_of_nodes() << std::endl;
        std::cout << "reduced_transformed_graph_edges=" << new_graph.number_of_edges() / 2 << std::endl;

        std::cout << "ils_time_to_best=" << ils_instance.get_last_update_time() << std::endl;
        std::cout << "total_time_to_best=" << ils_instance.get_last_total_update_time() << std::endl;
        std::cout << "clique_cover_size=" << near_linear_ils_cover << std::endl;
        std::cout << "verified_cover=no" << std::endl;
        std::cout << "optimal=" << (near_linear_ils_cover == partition_config.mis ? "yes":"unknown") << std::endl;
        }

        return 0;
    } else if (partition_config.run_type == "Redutransform") {

        timer total_timer;

        timer vcc_reductions_timer;
        redu_vcc oldVCC(G);
        std::vector<unsigned int> iso_degree;
        iso_degree.assign(G.number_of_nodes(), 0);
        std::vector<unsigned int> dom_degree;
        dom_degree.assign(G.number_of_nodes(), 0);
        reducer R;
        R.exhaustive_reductions(oldVCC, iso_degree, dom_degree);
        oldVCC.analyzeGraph(graph_filename, G, s, false /* don't check cover */);

        double vcc_reduction_time = vcc_reductions_timer.elapsed();

        oldVCC.buildKernel();

        graph_access kernel_graph;
        kernel_graph.start_construction(oldVCC.kernel_adj_list.size(), oldVCC.kernel_edges);
        for (NodeID v = 0; v < oldVCC.kernel_adj_list.size(); v++) {
            NodeID shadow_node = kernel_graph.new_node();
            kernel_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : oldVCC.kernel_adj_list[v]) {
                EdgeID shadow_edge = kernel_graph.new_edge(shadow_node, neighbor);
                kernel_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        kernel_graph.finish_construction();

        redu_vcc reduVCC(kernel_graph);

        timer transformation_timer;

        // Step 1: acyclic orientation
        std::vector<NodeID> order;
        order.reserve(G.number_of_nodes());
        for (NodeID node = 0; node < kernel_graph.number_of_nodes(); node++) {
            order.push_back(node);
        }

        // Step 2: line graph
        std::unordered_set<std::pair<NodeID,NodeID>> neighbor_hash;
        forall_nodes(kernel_graph, v) {
            forall_out_edges(kernel_graph, e, v) {
                NodeID u = kernel_graph.getEdgeTarget(e);
                neighbor_hash.insert(std::pair<NodeID,NodeID>(u,v));
                neighbor_hash.insert(std::pair<NodeID,NodeID>(v,u));
            } endfor
        } endfor

        std::unordered_map<std::pair<NodeID, NodeID>, NodeID> edge_to_vertex;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        std::vector<std::vector<NodeID>> edge_to_edges;

        std::vector<NodeID> const empty;

        std::cout << "Making filtered line graph..." << std::endl;
        forall_nodes(kernel_graph, v) {
            forall_out_edges(kernel_graph, e, v) {
                NodeID u = kernel_graph.getEdgeTarget(e);
                std::pair<NodeID,NodeID> node_pair(u, v);
                if (order[v] < order[u]) {
                    node_pair = std::pair<NodeID,NodeID>(v, u);
                }
                if (edge_to_vertex.find(node_pair) == edge_to_vertex.end()) {
                    vertex_to_edge.push_back(node_pair);
                    edge_to_vertex[node_pair] = vertex_to_edge.size() - 1;
                    edge_to_edges.push_back(empty);
                }
            } endfor

            forall_out_edges(kernel_graph, e1, v) {
                NodeID u = kernel_graph.getEdgeTarget(e1);
                std::pair<NodeID,NodeID> first_pair(v, u);
                if (order[u] < order[v]) 
                    first_pair = std::pair<NodeID,NodeID>(u, v);

                forall_out_edges(kernel_graph, e2, v) {
                    NodeID w = kernel_graph.getEdgeTarget(e2);
                    if (w <= u) continue;
                    std::pair<NodeID,NodeID> second_pair(v, w);
                    if (order[w] < order[v]) {
                        second_pair = std::pair<NodeID,NodeID>(w, v);
                    }

                    // Step 3: filter triples, TODO/DS: appears to work
                    if (std::get<0>(first_pair) == std::get<0>(second_pair) && 
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(first_pair),std::get<1>(second_pair))) != neighbor_hash.end()) && 
                        (neighbor_hash.find(std::pair<NodeID,NodeID>(std::get<1>(second_pair),std::get<1>(first_pair))) != neighbor_hash.end())) continue;

                    
                    edge_to_edges[edge_to_vertex[first_pair]].push_back(edge_to_vertex[second_pair]);
                    edge_to_edges[edge_to_vertex[second_pair]].push_back(edge_to_vertex[first_pair]);

                } endfor
            } endfor
        } endfor


        // Step 3: trim triples (integrated above)

        std::size_t cover_size_offset = R.get_cover_size_offset();
        std::cout << "cover_size_offset=" << cover_size_offset << std::endl;

        std::size_t num_edges = 0;
        for (std::vector<NodeID> & neighbors : edge_to_edges) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }

        double transformation_time = transformation_timer.elapsed();

        num_edges /= 2;

        graph_access new_graph;
        new_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        //std::cout << "Writing out transformed graph..." << std::endl;
        //graph_io::writeGraph(new_graph, "transformed.graph");

        // Step 4: Compute MIS
        if (true) // Exhaustive MIS Reductions + Sigmod + ILS
        {
        std::cout << "Making copy of transformed graph..." << std::endl;
        std::vector<std::vector<int>> line_copy;
        line_copy.reserve(edge_to_edges.size());
        for (std::vector<NodeID> const neighbors : edge_to_edges) {
            std::vector<int> int_neighbors(neighbors.size());
            for (std::size_t index = 0; index < neighbors.size(); index++) {
                int_neighbors[index] = neighbors[index];
            }
            line_copy.emplace_back(int_neighbors);
        }

        timer mis_reductions;

        std::cout << "Initializing branch-and-reduce..." << std::endl;
        branch_and_reduce_algorithm mis_bnr(line_copy, line_copy.size());
        std::cout << "Applying exhaustive MIS reductions..." << std::endl;
        mis_bnr.reduce();
        std::cout << "After exhaustive MIS reductions: " << line_copy.size() << " -> " << mis_bnr.number_of_nodes_remaining() << std::endl;

        std::vector<NodeID> reverse_mapping;

        double mis_reduction_time = mis_reductions.elapsed();

        graph_access full_kernel;
        mis_bnr.convert_adj_lists(full_kernel, reverse_mapping);
        std::cout << "Preparing ILS on full kernel of filtered line graph..." << std::endl;

        Graph mis_G;
        mis_G.read_graph(full_kernel);
        std::vector<std::vector<NodeID>> kernel;
        std::size_t offset = 0;
        std::cout << "Running reducing-peeling..." << endl;
        std::vector<bool> in_initial_is;

        timer peeling_timer;
        std::vector<NodeID> new_to_old_id_unused;
        unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id_unused, in_initial_is, offset);
        double peeling_time = peeling_timer.elapsed();
        std::cout << "Done with reducing-peeling..." << endl;
        assert(offset == 0);

        std::cout << "    Reduced " << full_kernel.number_of_nodes() << " -> " << kernel.size() << " nodes" << std::endl;
        std::cout << "    Offset =  " << offset << std::endl;

        std::size_t num_kernel_edges = 0;
        for (std::vector<NodeID> & neighbors : kernel) {
            std::sort(neighbors.begin(), neighbors.end());
            num_kernel_edges += neighbors.size();
        }
        num_kernel_edges /= 2;

        graph_access peeled_graph;
        peeled_graph.set_partition_count(2);
        peeled_graph.start_construction(kernel.size(), 2 * num_kernel_edges);
        for (NodeID v = 0; v < kernel.size(); v++) {
            NodeID shadow_node = peeled_graph.new_node();
            peeled_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
            peeled_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : kernel[v]) {
                EdgeID shadow_edge = peeled_graph.new_edge(shadow_node, neighbor);
                peeled_graph.setEdgeWeight(shadow_edge, 1);
            }
        }
        peeled_graph.finish_construction();

        //std::cout << "Writing out post-peeling graph peeled.graph" << std::endl;
        //graph_io::writeGraph(peeled_graph, "peeled.graph");

        std::size_t const ils_print_offset = kernel_graph.number_of_nodes() + cover_size_offset - mis_bnr.get_is_offset() - offset;

        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::cout << "Performing ILS..." << std::endl;
        ils_instance.perform_ils(ils_config, peeled_graph, ils_config.ils_iterations, total_timer.elapsed(), ils_print_offset, partition_config.mis);
        double time_for_ils = ils_timer.elapsed();
        cout << "ILS Running time: " << time_for_ils << std::endl;
        cout << "Total Running time (beginning to end): " << total_timer.elapsed() << std::endl;

        std::size_t const ils_cover = kernel_graph.number_of_nodes() - (ils_instance.get_best_solution_size() + mis_bnr.get_is_offset() + offset);

        std::cout << "Transform (Exhaustive+Peeling+ILS) Upper Bound = " << ils_cover + cover_size_offset << std::endl;

        std::size_t const clique_cover_size = ils_cover + cover_size_offset;

        std::cout << "BEGIN_OUTPUT_FOR_TABLES" << std::endl;
        std::cout << "run_type=" << partition_config.run_type << std::endl;
        std::cout << "vcc_reduction_time=" << vcc_reduction_time << std::endl;
        std::cout << "transformation_time=" << transformation_time << std::endl;
        std::cout << "transformed_graph_vertices=" << edge_to_edges.size() << std::endl;
        std::cout << "transformed_graph_edges=" << num_edges << std::endl;
        std::cout << "mis_reduction_time=" << mis_reduction_time << std::endl;
        std::cout << "reduced_transformed_graph_vertices=" << full_kernel.number_of_nodes() << std::endl;
        std::cout << "reduced_transformed_graph_edges=" << full_kernel.number_of_edges() / 2 << std::endl;
        std::cout << "peeling_time=" << peeling_time << std::endl;
        std::cout << "ils_time_to_best=" << ils_instance.get_last_update_time() << std::endl;
        std::cout << "total_time_to_best=" << ils_instance.get_last_total_update_time() << std::endl;
        std::cout << "clique_cover_size=" << clique_cover_size << std::endl;
        std::cout << "verified_cover=no" << std::endl;
        std::cout << "optimal=" << (clique_cover_size == partition_config.mis ? "yes":"unknown") << std::endl;
        }

        if (false) // Exhaustive MIS Reductions + ILS
        {
        std::cout << "Making copy of transformed graph..." << std::endl;
        std::vector<std::vector<int>> line_copy;
        line_copy.reserve(edge_to_edges.size());
        for (std::vector<NodeID> const neighbors : edge_to_edges) {
            std::vector<int> int_neighbors(neighbors.size());
            for (std::size_t index = 0; index < neighbors.size(); index++) {
                int_neighbors[index] = neighbors[index];
            }
            line_copy.emplace_back(int_neighbors);
        }

        std::cout << "Initializing branch-and-reduce..." << std::endl;
        branch_and_reduce_algorithm mis_bnr(line_copy, line_copy.size());
        std::cout << "Applying exhaustive MIS reductions..." << std::endl;
        mis_bnr.reduce();
        std::cout << "After exhaustive MIS reductions: " << line_copy.size() << " -> " << mis_bnr.number_of_nodes_remaining() << std::endl;

        std::vector<NodeID> reverse_mapping;
        graph_access full_kernel;
        mis_bnr.convert_adj_lists(full_kernel, reverse_mapping);
        std::cout << "Preparing ILS on full kernel of filtered line graph..." << std::endl;
        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = UINT_MAX;
        ils ils_instance;
        std::cout << "Performing ILS..." << std::endl;
        ils_instance.perform_ils(ils_config, full_kernel, ils_config.ils_iterations);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = kernel_graph.number_of_nodes() - (ils_instance.get_best_solution_size() + mis_bnr.get_is_offset());

        std::cout << "Transform (Full Kernel ILS) Upper Bound = " << ils_cover + cover_size_offset << std::endl;
        }

        if (false) // reducing-peeling initial IS only
        {
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        unsigned int res_mis = mis_G.degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction();

        std::size_t const cover_upper_bound = kernel_graph.number_of_nodes() - res_mis;
        std::cout << "Original Graph = " << kernel_graph.number_of_nodes() << ", MIS=" << res_mis << std::endl;
        std::cout << "Transform (Near Linear) Upper Bound = " << cover_upper_bound + cover_size_offset << std::endl;
        }

        if (false) // don't try any other variations
        {
        std::cout << "Building a sorted int version of adjacency list..." << std::endl;
        vector<vector<int>> adjlist(edge_to_edges.size());
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            std::sort(edge_to_edges[v].begin(), edge_to_edges[v].end());
            adjlist[v].reserve(edge_to_edges[v].size());
            for (NodeID const neighbor : edge_to_edges[v])
                adjlist[v].push_back(neighbor);
        }

        if (false) // Transform + ILS
        {
        std::size_t num_edges = 0;
        for (std::vector<NodeID> const &neighbors : edge_to_edges) {
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        std::cout << "Graph has " << edge_to_edges.size() << " vertices and " << num_edges << " edges" << std::endl;

        std::cout << "Constructing graph_access data structure..." << std::endl;

        graph_access new_graph;
        new_graph.start_construction(edge_to_edges.size(), 2 * num_edges);
        for (NodeID v = 0; v < edge_to_edges.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : edge_to_edges[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        std::cout << "Writing out transformed graph..." << std::endl;
        graph_io::writeGraph(new_graph, "transformed.graph");

        std::cout << "Preparing ILS..." << std::endl;
        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = ULONG_MAX;
        ils ils_instance;
        std::cout << "Performing ILS..." << std::endl;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = kernel_graph.number_of_nodes() - ils_instance.get_best_solution_size();

        std::cout << "Transform (ILS) Upper Bound = " << ils_cover + cover_size_offset << std::endl;
        }

        if (false) // BNR on transformed graph
        {
        branch_and_reduce_algorithm mis_bnr(adjlist, adjlist.size());

        timer bnr_timer;
        bool timeout = mis_bnr.solve(bnr_timer, partition_config.solver_time_limit) == -1;

        std::size_t const exact_cover = kernel_graph.number_of_nodes() - mis_bnr.get_current_is_size();

        std::cout << "BNR Running time: " << bnr_timer.elapsed() << std::endl;
        std::cout << "Transform (BNR) Upper Bound = " << exact_cover + cover_size_offset << std::endl;
        std::cout << "    Upper bound is " << (timeout ? "inexact" : "exact") << endl;
        }

        if (false) // reducing-peeling + ILS
        {
        Graph mis_G;
        mis_G.read_graph(edge_to_edges);
        std::vector<std::vector<NodeID>> kernel;
        std::size_t offset = 0;
        std::cout << "Computing Kernel..." << endl;
        std::vector<bool> in_initial_is;
        std::vector<NodeID> new_to_old_id_unused;
        unsigned int res_mis = mis_G.near_linear_kernel_and_offset(kernel, new_to_old_id_unused, in_initial_is, offset);
        std::cout << "    Reduced " << edge_to_edges.size() << " -> " << kernel.size() << " nodes" << std::endl;
        std::cout << "    Offset =  " << offset << std::endl;

        std::size_t num_edges = 0;
        for (std::vector<NodeID> & neighbors : kernel) {
            std::sort(neighbors.begin(), neighbors.end());
            num_edges += neighbors.size();
        }
        num_edges /= 2;

        graph_access new_graph;
        new_graph.set_partition_count(2);
        new_graph.start_construction(kernel.size(), 2 * num_edges);
        for (NodeID v = 0; v < kernel.size(); v++) {
            NodeID shadow_node = new_graph.new_node();
            new_graph.setPartitionIndex(shadow_node, in_initial_is[v] ? 1 : 0);
            new_graph.setNodeWeight(shadow_node, 1);
            for (NodeID const neighbor : kernel[v]) {
                EdgeID shadow_edge = new_graph.new_edge(shadow_node, neighbor);
                new_graph.setEdgeWeight(shadow_edge, 1);
            }
        }

        new_graph.finish_construction();

        timer kamis_timer;
        MISConfig kamis_config;
        kamis_config.seed = partition_config.seed;
        kamis_config.time_limit = partition_config.solver_time_limit;
        kamis_config.force_cand = 4;
        kamis_config.ils_iterations = 15000;

        timer ils_timer;
        MISConfig ils_config;
        ils_config.seed = partition_config.seed;
        ils_config.time_limit = partition_config.solver_time_limit;
        ils_config.force_cand = 4;
        ils_config.ils_iterations = ULONG_MAX;
        ils ils_instance;
        ils_instance.perform_ils(ils_config, new_graph, ils_config.ils_iterations);
        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const near_linear_ils_cover = kernel_graph.number_of_nodes() - (ils_instance.get_best_solution_size() + offset);

        std::cout << "Transform (Near Linear / ILS) Upper Bound = " << near_linear_ils_cover + cover_size_offset << std::endl;
        std::cout << "Done!" << std::endl;
        }
        } // other options
        return 0;
    } else if (partition_config.run_type == "test") {

        timer total_timer;
        timer transformation_timer;

        graph_access transformed_graph;
        std::vector<std::pair<NodeID, NodeID>> vertex_to_edge;
        transform_graph(G, transformed_graph, vertex_to_edge);

        double transformation_time = transformation_timer.elapsed();

        // Step 4: Compute MIS
        if (true) // just ILS
        {
        timer ils_timer;
        greedy_mis greedy;
        greedy.initial_partition(partition_config.seed, transformed_graph);

        ils ils_instance;
        run_ils(ils_instance, partition_config, transformed_graph, total_timer, G.number_of_nodes());

        cout << "ILS Running time: " << ils_timer.elapsed() << std::endl;

        std::size_t const ils_cover = G.number_of_nodes() - ils_instance.get_best_solution_size();

        std::cout << "Transform (ILS) Upper Bound = " << ils_cover << std::endl;

        std::vector<std::vector<NodeID>> clique_cover;
        size_t const num_cliques = G.number_of_nodes() - ils_instance.get_best_solution_size();

        transform_is_to_clique_cover(transformed_graph,
                                     vertex_to_edge,
                                     G.number_of_nodes(),
                                     ils_instance.get_best_solution(),
                                     num_cliques,
                                     clique_cover);
        }
        return 0;
    }

    redu_vcc reduVCC;
    branch_and_reduce B(G, reduVCC, partition_config);

    vertex_queue *queue = nullptr;
    if (partition_config.run_type == "cascading") queue = new vertex_queue(reduVCC);
    B.bandr(reduVCC, 0, queue, partition_config, s);
    B.analyzeGraph(graph_filename, G, reduVCC, s);

    // branch_and_reduce Bra(G, partition_config);
    // std::cout << "here" << std::endl;
    // std::vector<std::vector<NodeID>> cliques = Bra.sorted_enumerate(7, Bra.reduVCC.node_mis);
    // std::cout << std::endl;
    // for (std::vector<NodeID> &clique : cliques) {
    //   std::cout << "size: " << clique.size() << " mis: ";
    //   for (NodeID a : clique) {
    //     if (Bra.reduVCC.node_mis[a]) {
    //       std::cout << "1";
    //       break;
    //     }
    //   }
    //   std::cout << std::endl;
    // }

    // if (partition_config.run_type == "reductions_chalupa") {
    //   redu_vcc reduVCC(G);
    //   reducer R(G);
    //   // reduVCC.analyzeGraph(graph_filename, G, s);
    //   R.exhaustive_reductions(G, reduVCC);
    //   // reduVCC.analyzeGraph(graph_filename, G, s);
    //   reduVCC.build_cover(G);
    //   // reduVCC.analyzeGraph(graph_filename, G, s);
    //   reduVCC.solveKernel(G, partition_config, s);
    //   R.unwindReductions(G, reduVCC);
    //   reduVCC.analyzeGraph(graph_filename, G, s);
    //
    //   return 0;
    // }
    //
    // // branch_and_reduce B(G, partition_config);
    // if (partition_config.run_type == "brute") {
    //   branch_and_reduce B(G);
    //   B.brute_bandr(G, 0);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else if (partition_config.run_type == "reduMIS") {
    //   branch_and_reduce B(G, partition_config);
    //   B.reduMIS_bandr(G, 0);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else if (partition_config.run_type == "small_degree") {
    //   branch_and_reduce B(G, partition_config);
    //   B.small_degree_bandr(G, 0);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else if (partition_config.run_type == "sort_enum") {
    //   branch_and_reduce B(G, partition_config);
    //   B.sort_enum_bandr(G, 0, partition_config, s);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else if (partition_config.run_type == "chalupa_status") {
    //   branch_and_reduce B(G, partition_config);
    //   B.chalupa_status_bandr(G, 0, partition_config, s);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else if (partition_config.run_type == "KaMIS") {
    //   branch_and_reduce B(G);
    //   B.generate_mis_bandr(G, 0, partition_config, s);
    //   B.analyzeGraph(graph_filename, G, s);
    //   std::cout << "branches: " << B.branch_count << std::endl;
    // }
    // else {
    //   std::cout << "Error: required run-type" << std::endl;
    // }
    //
    // return 0;




}
