
#include <algorithm>
#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>

#include "redu_vcc.h"

void redu_vcc::getMIS(std::string file) {
  /* Generates node_mis mapping of minimum independent set from file */

  std::string line;

  std::ifstream mis_file (file);
  if (mis_file.is_open()) {
    while ( getline (mis_file, line)) {
      node_mis.push_back((int)line[0] - 48);
    }
    mis_file.close();
  }

  for (bool n : node_mis) if (n) { curr_mis++; };
}

redu_vcc::redu_vcc(graph_access &G, PartitionConfig &partition_config) : redu_structure(G) {

  if (!partition_config.mis_file.empty()) getMIS(partition_config.mis_file);
};


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

void redu_vcc::writeKernel(graph_access &G, std::string &filename) {
  if (remaining_nodes == 0) { return; }
  buildKernel(G);

  std::ofstream kernel_file;
  kernel_file.open (filename + "_kernel.txt");
  kernel_file << remaining_nodes << " " << kernel_edges / 2 << " 0\n";
  for (unsigned int i = 0; i < remaining_nodes; i++) {
    for (NodeID j : kernel_adj_list[i]) {
      kernel_file << j + 1 << " ";
    }
    kernel_file << "\n";
  }
  kernel_file.close();
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
