
#include <algorithm>
#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "redu_vcc.h"


void redu_vcc::solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t) {

  if (remaining_nodes == 0) { return; }

  buildKernel(G);

  cli *cli_instance;
  cli_instance = new cli(partition_config.seed, partition_config.mis);
  cli_instance->start_cli(kernel_adj_list, remaining_nodes, kernel_edges, t.elapsed(), partition_config.solver_time_limit);

  if (cli_instance->clique_cover.size() != 0){
      addKernelCliques(cli_instance->clique_cover);
  }
  else {
      std::cout << "Chalupa's algorithm unable to solve in given time." << std::endl;
  }

  delete(cli_instance);

}

void redu_vcc::analyzeGraph(std::string &filename, graph_access &G, timer &t){

    std::cout << filename << ", ";

    std::cout << G.number_of_nodes() << ", ";
    std::cout << G.number_of_edges() << ", ";

    std::cout <<  remaining_nodes << ", ";

    std::cout << t.elapsed() << ", ";

    std::cout << clique_cover.size() << std::endl;

    validateCover(G);

}
