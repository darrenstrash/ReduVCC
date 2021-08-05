
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

#include "mis/kernel/branch_and_reduce_algorithm.h"

#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_vcc/redu_vcc.h"
#include "redu_vcc/reducer.h"
#include "branch_and_reduce/b_and_r.h"

#include "sigmod_mis/Graph.h"
//
// #include "mis/ils/ils.h"
// #include "mis/mis_config.h"

// #include "redu_vcc/graph_vcc.h"

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


    redu_vcc reduVCC;
    branch_and_reduce B(G, reduVCC, partition_config);

    vertex_queue *queue = nullptr;
    if (partition_config.run_type == "edge") { B.edge_bandr(reduVCC, 0, queue, partition_config, s, 0); }
    else { B.bandr(reduVCC, 0, queue, partition_config, s);}
    B.analyzeGraph(graph_filename, G, reduVCC, s);

    return;

    // if (partition_config.run_type == "Redu") {
    //     redu_vcc reduVCC(G);
    //     std::vector<unsigned int> iso_degree;
    //     iso_degree.assign(G.number_of_nodes(), 0);
    //     std::vector<unsigned int> dom_degree;
    //     dom_degree.assign(G.number_of_nodes(), 0);
    //     reducer R;
    //     R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);
    //     reduVCC.analyzeGraph(graph_filename, G, s);
    //     return 0;
    // }
    // if (partition_config.run_type == "ReduVCC") {
    //     redu_vcc reduVCC(G);
    //     std::vector<unsigned int> iso_degree;
    //     iso_degree.assign(G.number_of_nodes(), 0);
    //     std::vector<unsigned int> dom_degree;
    //     dom_degree.assign(G.number_of_nodes(), 0);
    //     reducer R;
    //     R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);
    //     reduVCC.analyzeGraph(graph_filename, G, s);
    //     reduVCC.build_cover();
    //     reduVCC.solveKernel(partition_config, s);
    //     R.unwindReductions(reduVCC);
    //     reduVCC.analyzeGraph(graph_filename, G, s);
    //     return 0;
    // }
    //
    // redu_vcc reduVCC;
    // branch_and_reduce B(G, reduVCC, partition_config);
    //
    // vertex_queue *queue = nullptr;
    // if (partition_config.run_type == "cascading") queue = new vertex_queue(reduVCC);
    // B.bandr(reduVCC, 0, queue, partition_config, s);
    // B.analyzeGraph(graph_filename, G, reduVCC, s);





}
