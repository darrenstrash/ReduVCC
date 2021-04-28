
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

void redu_vcc::generateAdjList(graph_access &G) {
  /* Generates adjacency list from graph */

  forall_nodes(G, v){
      std::vector<NodeID> N_v;    // open neighborhood of v

      forall_out_edges(G, e, v){
          NodeID u = G.getEdgeTarget(e);
          N_v.push_back(u);

      } endfor
      adj_list.push_back(N_v);

  } endfor
}

void redu_vcc::subgraph_map(std::vector<NodeID> &subgraph_nodes) {
  // maps subgraph vertices from parent to new instance

  parent_to_self_map.resize(num_nodes);
  self_to_parent_map.resize(num_nodes);
  num_nodes = 0;

  for (NodeID old_v : subgraph_nodes) {
    parent_to_self_map[old_v] = num_nodes;
    self_to_parent_map[num_nodes] = old_v;

    num_nodes++;
  }
}

void redu_vcc::generateAdjList(redu_vcc* parent) {
  // generates adj_list from parent

  for (NodeID v = 0; v < num_nodes; v++){
    // get v from parent
    NodeID old_v = self_to_parent_map[v];
    // std::cout << old_v << std::endl;

    std::vector<NodeID> adj;
    for (NodeID old_u : parent->adj_list[old_v]) {
      if (!parent->node_status[old_u]) continue;
      NodeID u = parent_to_self_map[old_u];

      adj.push_back(u);
      // std::cout << u;
    }
    std::sort(adj.begin(), adj.end());

    adj_list.push_back(adj);
  }
}

redu_vcc::redu_vcc(graph_access &G, PartitionConfig &partition_config) {

  num_nodes = G.number_of_nodes();
  generateAdjList(G);

  init();

  if (!partition_config.mis_file.empty()) getMIS(partition_config.mis_file);
};

redu_vcc::redu_vcc(redu_vcc* parent, std::vector<NodeID> &subgraph_nodes) {

  parent_to_self_map.resize(parent->num_nodes);
  self_to_parent_map.resize(parent->num_nodes);

  subgraph_map(subgraph_nodes);
  generateAdjList(parent);

  init();
  if (!partition_config.mis_file.empty()) getMIS(partition_config.mis_file);

}

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
