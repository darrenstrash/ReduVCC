
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

#include "redu_vcc/reducer.h"
#include "branch_and_reduce/b_and_r.h"

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


    branch_and_reduce B(G);
    B.getMIS(partition_config.mis_file);
    // B.branch(G, 0);
    // B.prune_branch(G, 0, B.mis);
    B.small_deg_branch(G, 0, B.mis);
    B.analyzeGraph(graph_filename, G, t);

    // // std::cout << "io time: " << t.elapsed()  << std::endl;
    //
    // timer s;
    // // B.enumerate(0);
    // if (partition_config.run_type == "brute"){
    //   B.brute(G);
    // }
    // else if (partition_config.run_type == "prune") {
    //   B.prune(G);
    // }
    // else if (partition_config.run_type == "small_degree") {
    //   B.min_degree_prune(G);
    // }
    // else { return; }
    // // unsigned int i = 0;
    // // B.exhaustive_reductions(G, i, i);
    // B.analyzeGraph(graph_filename, G, s);

     // reducer R;
     // R.init(G);
     // std::vector<unsigned int> iso_stats = R.bruteISO(G);
     // R.analyzeGraph(graph_filename, G, s);
     // std::vector<unsigned int> d2_stats = R.bruteD2(G);
     // R.analyzeGraph(graph_filename, G, s);
     // std::vector<unsigned int> twin_stats = R.bruteTWIN(G);
     // R.analyzeGraph(graph_filename, G, s);
     // std::vector<unsigned int> dom_stats = R.bruteDOM(G);
     // R.analyzeGraph(graph_filename, G, s);
     // std::vector<unsigned int> crown_stats = R.bruteCROWN(G);
     // R.analyzeGraph(graph_filename, G, s);
     // R.analyzeGraph(graph_filename, G, s);
     // R.buildCover(G);
     // R.solveKernel(G, partition_config, s);
     // R.unwindReductions(G);
     // R.analyzeGraph(graph_filename, G, s);
     // R.undoReductions(G, l);
     // // // // // //
     // R.analyzeGraph(graph_filename, G, s);

}
