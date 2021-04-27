
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


std::vector<std::vector<NodeID>> redu_vcc::generateAdjList(graph_access &G) {
  // Generates adjacency list from graph

  std::vector<std::vector<NodeID>> adj_list;

  forall_nodes(G, v){
      std::vector<NodeID> N_v;    // open neighborhood of v

      forall_out_edges(G, e, v){
          NodeID u = G.getEdgeTarget(e);
          N_v.push_back(u);

      } endfor
      adj_list.push_back(N_v);

  } endfor

  return adj_list;
}

unsigned int redu_vcc::subgraph_map(std::vector<NodeID> &subgraph_nodes) {
  // maps subgraph vertices from parent to new instance

  unsigned int num_nodes = 0;

  for (NodeID old_v : subgraph_nodes) {
    parent_to_self_map[old_v] = num_nodes;
    self_to_parent_map[num_nodes] = old_v;

    num_nodes++;
  }

  return num_nodes;
}

std::vector<std::vector<NodeID>> redu_vcc::generateAdjList(redu_vcc* parent, unsigned int &num_nodes) {
  // generates adj_list from parent

  std::vector<std::vector<NodeID>> adj_list;

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

    adj_list.push_back(adj);
  }

  return adj_list;
}


redu_vcc::redu_vcc(graph_access &G) {

  unsigned int num_nodes = G.number_of_nodes();
  std::vector<std::vector<NodeID>> adj_list = generateAdjList(G);

  init(num_nodes, adj_list);

}

redu_vcc::redu_vcc(redu_vcc* parent, std::vector<NodeID> &subgraph_nodes) {

  parent_to_self_map.resize(parent->num_nodes);
  self_to_parent_map.resize(parent->num_nodes);

  unsigned int num_nodes = subgraph_map(subgraph_nodes);
  std::vector<std::vector<NodeID>> adj_list = generateAdjList(parent, num_nodes);

  init(num_nodes, adj_list);

}

std::vector<NodeID> redu_vcc::find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining) {

  std::vector<NodeID> current_nodes;
  std::vector<NodeID> queue;

  NodeID v = 0;
  while (visited_nodes[v]) v ++;
  visited_nodes[v] = true; visit_remaining--;
  current_nodes.push_back(v);

  queue.push_back(v);

  while (!queue.size() == 0) {

    v = queue.front();
    queue.erase(queue.begin());


    for (NodeID u : adj_list[v]) {
      if (!visited_nodes[u]) {
        visited_nodes[u] = true; visit_remaining--;
        current_nodes.push_back(u);
        queue.push_back(u);
      }
    }

  }

  // for (NodeID a : current_nodes) {
  //   std::cout << a << ", ";
  // }
  // std::cout << std::endl;

  return current_nodes;
}

void redu_vcc::decompose_components() {

  std::vector<bool> visited_nodes;
  for (bool status : node_status) visited_nodes.push_back(!status);
  unsigned int visit_remaining = node_status.size();

  std::vector<NodeID> subgraph_nodes = find_component(visited_nodes, visit_remaining);
  if (visit_remaining == 0) {
    return;
  }


  redu_vcc child(this, subgraph_nodes);
  children.push_back(child);

  while (visit_remaining > 0) {
    subgraph_nodes = find_component(visited_nodes, visit_remaining);
    child = redu_vcc(this, subgraph_nodes);
    children.push_back(child);

  }
}

void redu_vcc::merge_child(redu_vcc& child) {

}

void redu_vcc::merge_kernel() {

}

void redu_vcc::build_cover() {

  // builds existing cliques
  if (remaining_nodes != 0){
    if (children.size() == 0) {
      merge_kernel();
    }
    else {
      for (redu_vcc& child : children) {
        // creates complete cover of child
        child.build_cover();
        // integrates child cover
        merge_child(child);
      }
    }
  }

  for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(this);


}


// void redu_vcc::getMIS(std::string file) {
//   /* Generates node_mis mapping of minimum independent set from file */
//
//   std::string line;
//
//   std::ifstream mis_file (file);
//   if (mis_file.is_open()) {
//     while ( getline (mis_file, line)) {
//       node_mis.push_back((int)line[0] - 48);
//     }
//     mis_file.close();
//   }
//
//   for (bool n : node_mis) if (n) { curr_mis++; };
// }
//
// redu_vcc::redu_vcc(graph_access &G, PartitionConfig &partition_config) : redu_structure(G) {
//
//   if (!partition_config.mis_file.empty()) getMIS(partition_config.mis_file);
// };
//
//
//
//
//
//
// void redu_vcc::solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t) {
//
//   if (remaining_nodes == 0) { return; }
//
//   buildKernel(G);
//
//   cli *cli_instance;
//   cli_instance = new cli(partition_config.seed, partition_config.mis);
//   cli_instance->start_cli(kernel_adj_list, remaining_nodes, kernel_edges, t.elapsed(), partition_config.solver_time_limit);
//
//   if (cli_instance->clique_cover.size() != 0){
//       addKernelCliques(cli_instance->clique_cover);
//   }
//   else {
//       std::cout << "Chalupa's algorithm unable to solve in given time." << std::endl;
//   }
//
//   delete(cli_instance);
//
// }
//
// void redu_vcc::analyzeGraph(std::string &filename, graph_access &G, timer &t){
//
//     std::cout << filename << ", ";
//
//     std::cout << G.number_of_nodes() << ", ";
//     std::cout << G.number_of_edges() << ", ";
//
//     std::cout <<  remaining_nodes << ", ";
//
//     std::cout << t.elapsed() << ", ";
//
//     std::cout << clique_cover.size() << std::endl;
//
//     validateCover(G);
//
// }
//
//
// redu_vcc redu_vcc::subgraph(std::vector<NodeID> subgraph_nodes) {
//
//     /* Needs:
//         - create a back and forth mapping
//         -
//     */
// }
