
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

#include "mis/ils/ils.h"
#include "mis/mis_config.h"

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

    MISConfig config;
    config.console_log = true;
    config.time_limit = 60;
    config.force_cand = 10;
    ils new_ils;
    new_ils.perform_ils(config, G, 10);	
   return 0;

	
    if (partition_config.run_type == "bandr") {
      branch_and_reduce B(G, partition_config);
      // B.getMIS(partition_config.mis_file);
      // B.branch(G, 0);
      // B.prune_branch(G, 0, B.mis);
      // B.small_deg_branch(G, 0, B.mis);
      // B.lower_bound_branch(G, 0);
      B.sort_enumerate_branch(G, 0);
      B.analyzeGraph(graph_filename, G, s);
      return 0;
    }

    redu_vcc reduVCC(G);
    reducer R(G);
    // reduVCC.analyzeGraph(graph_filename, G, s);
    R.exhaustive_reductions(G, reduVCC);
    // reduVCC.analyzeGraph(graph_filename, G, s);
    reduVCC.build_cover(G);
    // reduVCC.analyzeGraph(graph_filename, G, s);
    reduVCC.solveKernel(G, partition_config, s);
    R.unwindReductions(G, reduVCC);
    reduVCC.analyzeGraph(graph_filename, G, s);

    return 0;




}
