
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

void redu_vcc::generateAdjList(redu_vcc *parent, std::vector<NodeID> &child_nodes) {
  /* Generates adjacency list from parent graph_vcc instance */

  // assign nodes in graph
  num_nodes = child_nodes.size();

  // std::cout << "begin map" << std::endl;
  // mapping from parent to child
  std::vector<NodeID> parent_to_self;
  std::vector<bool> subgraph_node; // marks vertices that are subgraph nodes
  subgraph_node.assign(parent->num_nodes, false);
  parent_to_self.resize(parent->num_nodes);
  self_to_parent.resize(child_nodes.size()); // allocate to size of original parent
  NodeID v = 0;

  for (NodeID u : child_nodes) {
    parent_to_self[u] = v;
    subgraph_node[u] = true;
    self_to_parent[v] = u;

    v++;
  }

  // std::cout << "begin trans" << std::endl;

  for (NodeID u : child_nodes) {
    std::vector<NodeID> N_v;
    for (NodeID w : parent->adj_list[u]) {
      if (!subgraph_node[w]) continue;
      NodeID child_w = parent_to_self[w];
      N_v.push_back(child_w);
    }
    std::sort(N_v.begin(), N_v.end());
    adj_list.push_back(N_v);
  }

    // std::cout << "end" << std::endl;
}

void redu_vcc::init() {

  // assign status of nodes
  node_status.assign(num_nodes, true);
  fold_node.assign(num_nodes, false);
  remaining_nodes = num_nodes;
  // allocate for graph cover
  node_clique.resize(num_nodes);

  // initialize mis mapping to 0
  curr_mis = 0;

  // allocate two scratch vectors
  scratch1.assign(num_nodes, false);
  scratch2.assign(num_nodes, false);

  // assign first cliqueID to 0
  next_cliqueID = 0;
}

redu_vcc::redu_vcc(graph_access &G) {

  num_nodes = G.number_of_nodes();
  // produce adjacency list
  generateAdjList(G);

  init();

}

redu_vcc::redu_vcc(redu_vcc *parent, std::vector<NodeID> &child_nodes) {

  num_nodes = child_nodes.size();
  generateAdjList(parent, child_nodes);

  init();
}

std::vector<NodeID> redu_vcc::find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining) {

  std::vector<NodeID> current_nodes;
  std::vector<NodeID> queue;
  unsigned int queue_size = 0;

  NodeID v = 0;
  while (visited_nodes[v]) v ++;
  visited_nodes[v] = true; visit_remaining--;
  current_nodes.push_back(v);

  queue.push_back(v);
  // if (v == 1138517) std::cout << "HERE" << std::endl;
  queue_size++;

  // std::cout << "ready for queue" << std::endl;

  while (!queue.size() == 0) {
    // std::cout << queue.size() << ", " << queue_size << std::endl;
    v = queue.front();
    queue.erase(queue.begin());
    queue_size--;

    // std::cout << "queue front: " << v << std::endl;
    // std::cout << " adj_size: " << adj_list[v].size() << std::endl;
    // std::cout << "added to queue: ";
    for (NodeID u : adj_list[v]) {
      if (!visited_nodes[u]) {
        visited_nodes[u] = true; visit_remaining--;
        current_nodes.push_back(u);
        // std::cout << u << ", ";
        // if (u == 1138517) std::cout << "HERE" << std::endl;
        queue.push_back(u);
        queue_size++;
      }
    }
    // std::cout << std::endl;
  }

  // std::cout << "comp found" << std::endl;

  // for (NodeID a : current_nodes) {
  //   std::cout << a << ", ";
  // }
  // std::cout << std::endl;

  return current_nodes;
}

std::vector<redu_vcc> redu_vcc::decompose() {

  // printAdjList();

  std::vector<redu_vcc> children;

  std::vector<bool> visited_nodes;
  for (bool status : node_status) visited_nodes.push_back(!status);
  unsigned int visit_remaining = remaining_nodes;

  std::vector<NodeID> subgraph_nodes = find_component(visited_nodes, visit_remaining);
  // std::cout << "num subgraph: "<< subgraph_nodes.size() << std::endl;
  if (visit_remaining == 0) {
    // std::cout << "no comp" << std::endl;
    return children;
  }

  // for (NodeID a : subgraph_nodes) {
  //   std::cout << a << ", ";
  // }
  // std::cout << std::endl;

  redu_vcc child(this, subgraph_nodes);
  children.push_back(child);
  // std::cout << "child created" << std::endl;

  while (visit_remaining > 0) {
    subgraph_nodes = find_component(visited_nodes, visit_remaining);
    // for (NodeID a : subgraph_nodes) {
    //   std::cout << a << ", ";
    // }
    // std::cout << std::endl;
    // std::cout << subgraph_nodes.size() << std::endl;
    child = redu_vcc(this, subgraph_nodes);
    children.push_back(child);
  }

  // std::cout << "all children" << std::endl;

  return children;

}

void redu_vcc::addCliquesToParent(redu_vcc &parent) {

  for (std::vector<NodeID> &clique : clique_cover) {
    std::vector<NodeID> parent_clique;

    for (NodeID &v : clique) {
      NodeID old_v = self_to_parent[v];
      parent_clique.push_back(old_v);
    }

    std::sort(parent_clique.begin(), parent_clique.end());
    parent.addCliqueToCover(parent_clique);
  }
}


void redu_vcc::build_cover(){
  /* Constructs clique cover from node_clique mapping. */

  // clears clique_cover
  clique_cover.clear();
  clique_cover.assign(next_cliqueID, {});

  for (NodeID v = 0; v < num_nodes; v++) {
  // forall_nodes(G, v) {
    if (fold_node[v]) { continue; } // if node in fold, skip
    if (node_status[v]) { continue; } // if node still in graph, skip

    unsigned int cliqueID = node_clique[v];
    clique_cover[cliqueID].push_back(v);
  }
  // } endfor

  // prepare to solve, by setting solve node_clique mapping and next cliqueID
  solve_node_clique = node_clique;
  next_solvecliqueID = next_cliqueID;
}

bool redu_vcc::cliqueInG(graph_access &G, std::vector<NodeID> &clique) {
  /* Tests if clique is valid in G */

  for (unsigned int i = 0; i < clique.size() -1; i++) {
    NodeID v = clique[i];
    unsigned int j = i + 1;
    forall_out_edges(G, e, v){
        NodeID u = G.getEdgeTarget(e);

        if (j >= clique.size()) { break; }
        else if (u == clique[j]) {
          j++;
          continue;
        }
        else if (u > clique[j]) { return false;}

    } endfor
  }

  return true;
}

void redu_vcc::validateCover(graph_access &G) {

  std::vector<bool> temp_status;
  temp_status.assign(G.number_of_nodes(), true);

  for (std::vector<NodeID> clique : clique_cover ) {
    if (clique.size() == 0) { std::cout << "Null clique" << std::endl; return; }

    for (NodeID v : clique) {
      // std::cout << v << ", ";
      if (temp_status[v] == false) {
        std::cout << "Overlap" << std::endl;
        return;
      }
      else { temp_status[v] = false; }
    }

    // std::cout << std::endl;

    if (!cliqueInG(G, clique)) {
      printVectorSet(clique);
      std::cout << "Invalid clique" << std::endl;
      return; }
    //
    // std::cout << std::endl;
  }
}

void redu_vcc::assignMaps() {

  old_to_new_map.clear();
  old_to_new_map.resize(num_nodes);
  new_to_old_map.clear();
  new_to_old_map.resize(num_nodes);

  int j = 0;
  for (unsigned int i = 0; i < num_nodes; i++){
      if (!node_status[i]){
          continue;
      }
      old_to_new_map[i] = j;
      new_to_old_map[j] = i;
      j++;
  }
}

void redu_vcc::buildKernel() {

  assignMaps();

  kernel_adj_list.clear();
  kernel_adj_list.resize(remaining_nodes);
  kernel_edges = 0;

  for (unsigned int i = 0; i < num_nodes; i++) {
    if (!node_status[i]) { continue; }

    int new_v = old_to_new_map[i];

    std::vector<int> adj;
    for (unsigned int j = 0; j < adj_list[i].size(); j++){
        NodeID u  = adj_list[i][j];
        if (!node_status[u]) { continue; }

        int new_u = old_to_new_map[u];
        adj.push_back(new_u);
        kernel_edges++;
    }
    std::sort(adj.begin(), adj.end());
    kernel_adj_list[new_v] = adj;
  }
}

void redu_vcc::addKernelCliques(std::vector<std::vector<int>> &clique_set){

  for (unsigned int i = 0; i < clique_set.size(); i++){
      std::vector<NodeID> clique;

      for (unsigned int j = 0; j < clique_set[i].size(); j++){
          int v = clique_set[i][j];
          NodeID old_v = new_to_old_map[v];

          solve_node_clique[old_v] = false;
          clique.push_back(old_v);
      }
      std::sort(clique.begin(), clique.end());

      addCliqueToCover(clique);
  }
}

void redu_vcc::addCrownCliques(std::vector<std::vector<NodeID>> &crown_cliques, std::vector<std::vector<int>> &clique_set) {

  for (unsigned int i = 0; i < clique_set.size(); i++){
      std::vector<NodeID> clique;

      for (unsigned int j = 0; j < clique_set[i].size(); j++){
          int v = clique_set[i][j];
          NodeID old_v = new_to_old_map[v];
          // std::cout << old_v << std::endl;
          // solve_node_clique[old_v] = false;
          clique.push_back(old_v);
      }
      std::sort(clique.begin(), clique.end());

      // printVectorSet(clique);
      addClique(clique);
      removeVertexSet(clique);

      crown_cliques.push_back(clique);
  }

}

unsigned int redu_vcc::adj_size(NodeID v) {

  unsigned int size = 0;
  for (NodeID u: adj_list[v]) {
    if (node_status[u]) { size++; }
  }

  return size;
}

std::vector<NodeID> redu_vcc::curr_adj_list(NodeID v) {

  std::vector<NodeID> curr_adj_list;

  for (NodeID u : adj_list[v]){
    if (node_status[u]) { curr_adj_list.push_back(u); }
  }

  return curr_adj_list;
}


void redu_vcc::removeVertex(NodeID v) {
  // removes vertex

  node_status[v] = false;
  remaining_nodes--;

  if (!node_mis.empty() && node_mis[v]) curr_mis--;
}

void redu_vcc::addVertex(NodeID v) {
  // adds vertex

  node_status[v] = true;
  remaining_nodes++;

  if (!node_mis.empty() && node_mis[v]) curr_mis++;
}

// void redu_vcc::addClique(std::vector<NodeID> &clique) {
//
//   std::sort(clique.begin(), clique.end());
//   clique_cover.push_back(clique);
//   unsigned int cliqueID = clique_cover.size() - 1;
//
//   for (NodeID u : clique){
//       node_clique[u] = cliqueID;
//   }
//
// }

void redu_vcc::addClique(std::vector<NodeID> &clique) {

  for (NodeID u : clique) {
    node_clique[u] = next_cliqueID;
  }

  next_cliqueID++;
}

void redu_vcc::addCliqueToCover(std::vector<NodeID> &clique) {
  // adds clique so solve node structure

  for (NodeID u : clique) {
    solve_node_clique[u] = next_solvecliqueID;
  }
  clique_cover.push_back(clique);
  next_solvecliqueID++;
  // std::cout << next_solvecliqueID << std::endl;
}



// std::vector<NodeID> redu_vcc::pop_clique() {
//
//   std::vector<NodeID> clique = clique_cover.back();
//   clique_cover.pop_back();
//   return clique;
// }

void redu_vcc::pop_clique(std::vector<NodeID> &clique) {

  for (NodeID u : clique) {
    node_clique[u] = node_clique.size();
  }

  next_cliqueID--;
}

unsigned int redu_vcc::getCliqueID(NodeID &v) {

  return node_clique[v];
}

std::vector<NodeID> redu_vcc::getClique(NodeID &v) {

  unsigned int cliqueID = getCliqueID(v);
  std::vector<NodeID> clique = clique_cover[cliqueID];

  return clique;
}

void redu_vcc::replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique){

  std::sort(new_clique.begin(), new_clique.end());

  for (NodeID a : new_clique) {
    solve_node_clique[a] = cliqueID;
  }
  clique_cover[cliqueID] = new_clique;
}


void redu_vcc::removeVertexSet(std::vector<NodeID> &S) {

  for (NodeID v : S) { removeVertex(v); };
}

void redu_vcc::addVertexSet(std::vector<NodeID> &S) {
  for (NodeID v : S) { addVertex(v); };
}

// void redu_vcc::clearScratch(std::vector<bool> &scratch) {
//
//   std::fill(scratch.begin(), scratch.end(), false);
// }

void redu_vcc::printAdjList() {

  for (unsigned int i = 0; i < node_status.size(); i++) {
    if (!node_status[i]) { continue; }
    std::cout << "N(" << i << "): [";

    for (NodeID u : adj_list[i]) {
      if (!node_status[u]) { continue; }
      std::cout << u << ", ";
    }
    std::cout << "]" << std::endl;
  }
}

void redu_vcc::printAdjList(NodeID v) {

  if (!node_status[v]) {
    std::cout << v << " /notin adj_list" << std::endl;
    return;
  }

  std::cout << "N(" << v << "): [";
  for (NodeID u : adj_list[v]) {
    if (!node_status[u]) { continue; }
    std::cout << u << ", ";
  }
  std::cout << "]" << std::endl;
}

void redu_vcc::printNeighborhood(NodeID v) {

  printAdjList(v);
  for (NodeID u : adj_list[v]) { std::cout << "  "; printAdjList(u); }
}

void redu_vcc::printVectorSet(std::vector<NodeID> S){

  for (unsigned int i = 0; i < S.size(); i++){
    std::cout << S[i] << ", ";
  }
  std::cout << std::endl;
}



//
// #include <algorithm>
// #include <argtable3.h>
// #include <iostream>
// #include <fstream>
// #include <math.h>
// #include <regex.h>
// #include <sstream>
// #include <stdio.h>
// #include <string.h>
//
// #include "redu_vcc.h"
//
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

redu_vcc::redu_vcc(graph_access &G, PartitionConfig &partition_config) : redu_vcc(G) {

  if (!partition_config.mis_file.empty()) getMIS(partition_config.mis_file);
};


void redu_vcc::solveKernel(PartitionConfig &partition_config, timer &t) {

  if (remaining_nodes == 0) { return; }

  buildKernel();

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
