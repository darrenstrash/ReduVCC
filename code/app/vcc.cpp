
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
     std::cout << "io time: " << t.elapsed()  << std::endl;

    timer s;

     reducer R(G);
     R.bruteISO(G);
     R.bruteD2(G);
     R.bruteTWIN(G);
     R.bruteDOM(G);
     R.bruteCROWN(G);
     R.analyzeGraph(graph_filename, G, s);
     R.solveKernel(G, partition_config, s);
     R.unwindReductions(G);

     R.analyzeGraph(graph_filename, G, s);

   //  std::cout << "here" << std::endl;
   //
   //  timer s;
   //
   //  scratch1.assign(G.number_of_nodes(), false);
   //  scratch2.assign(G.number_of_nodes(), false);
   //
   //  Reducer R(G);
   //
   //  if (partition_config.run_type == "reduction"){
   //    R.performReductions(G);
   //    R.analyzeGraph(graph_filename, G, s);
   //    std::cout << std::endl;
   //  }
   //
   // R.solveKernel(G, partition_config, s);
   // R.unwindReductions(G);
   // R.analyzeGraph(graph_filename, G, s);
   // R.validCover(G);

}
