
#include <algorithm>
#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "redu_structure.h"


void redu_structure::init(unsigned int node_count, std::vector<std::vector<NodeID>> &adjacency) {
  // initializes attributs of redu_structure

  num_nodes = node_count;
  adj_list = adjacency;

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

//
// void redu_structure::build_child_cover() {
//
//   for (redu_vcc &child : children) {
//     // build complete cover of child
//     child.build_cover();
//
//     // map child cliques to parent
//     for (std::vector<NodeID> &child_clique : child.clique_cover) {
//       std::vector<NodeID> clique;
//       for (NodeID &child_v : child_clique) {
//         NodeID v = child.self_to_parent_map[child_v];
//         clique.push_back(v);
//         node_clique[v] = next_cliqueID;
//       }
//       addCliqueToCover(clique);
//     }
//
//   }
// }
//
// void redu_structure::kernelAdjList(unsigned int &kernel_edges) {
//   // build new int adj_list
//
//   std::vector<std::vector<int> kernel_adj_list;
//
//   for (NodeID v = 0; v < num_nodes; v++) {
//     int new_v = self_to_chalupa_map[v];
//     std::vector<int> new_adj;
//     for (NodeID u : adj_list[v]) {
//       NodeID new_u = self_to_chalupa_map[u];
//       new_adj.push_back(new_u);
//       kernel_edges++;
//     }
//     std::sort(new_adj.begin(), new_adj.end());
//     kernel_adj_list.push_back(new_adj)
//   }
//
//   return kernel_adj_list
// }
//
// void redu_structure::solveKernel() {
//   // solve kernel using Chalupa's IG algorithm
//
//   // allocate space for maps
//   std::vector<int> self_to_chalupa_map.resize(num_nodes);
//   std::vector<NodeID chalupa_to_self_map.resize(remaining_nodes);
//
//   // assign maps
//   int new_v = 0;
//   for (NodeID v = 0; v < num_nodes; v++) {
//     if (!node_status[v]) continue;
//     self_to_chalupa_map[v] = new_v;
//     chalupa_to_self_map[new_v] = v;
//   }
//
//   unsigned int kernel_edges;
//   std::vector<std::vector<int> kernel_adj_list = kernelAdjList(kernel_edges);
//
//   cli *cli_instance;
//   cli_instance = new cli(partition_config.seed, partition_config.mis);
//   cli_instance->start_cli(kernel_adj_list, remaining_nodes, kernel_edges, t.elapsed(), partition_config.solver_time_limit);
//
//   if (cli_instance->clique_cover.size() != 0) addKernelCliques(cli_instance->clique_cover);
//   else std::cout << "Chalupa's algorithm unable to solve in given time." << std::endl;
//   delete(cli_instance);
//
//
// }
//
// void redu_structure::build_cover(){
//   // Constructs clique cover from node_clique mapping.
//
//   // clears clique_cover
//   clique_cover.clear();
//   clique_cover.assign(next_cliqueID, {});
//
//
//   for (NodeID v = 0; v < num_nodes; v++)
//     if (fold_node[v]) { continue; } // if node in fold, skip
//     if (node_status[v]) { continue; } // if node still in graph, skip
//
//     unsigned int cliqueID = node_clique[v];
//     clique_cover[cliqueID].push_back(v);
//   }
//
//   // prepare to solve, by setting solve node_clique mapping and next cliqueID
//   solve_node_clique = node_clique;
//   next_solvecliqueID = next_cliqueID;
//
//   if (remaining_nodes > 0) {
//     if (children.size() != 0) build_child_cover();
//     else chalupa;
//   }
//
//   // unwind reductions
// }
//
// bool redu_structure::safeEdge(NodeID &v, NodeID &u) {
//   // checks if edge is valid in original instance
//
//   for (NodeID &w : adj_list[v]) {
//     if (w == u) return true;
//   }
//   return false;
// }
//
// bool redu_structure::safeClique(std::vector<NodeID> &clique) {
//   // checks to see if a clique is in the unreduced instance
//
//   for (unsigned int i = 0; i < clique.size() - 1; i++) {
//     NodeID v = clique[i];
//
//     for (unsigned int j = i + 1; j < clique.size(); j++) {
//       NodeID u = clique[j];
//       if (!safeEdge(v, u)) return false;
//     }
//   }
//
//   return true;
// }
//
// void redu_structure::validateCover() {
//   // tests if current cover is vaild for original instance
//
//   unsigned int check_nodes = num_nodes;
//   std::vector<bool> check_status;
//   check_status.assign(num_nodes, true);
//
//   for (std::vector<NodeID> &clique: clique_cover) {
//     if (clique.size() == 0) {
//       std::cout << "Null clique" << std::endl;
//       return;
//     }
//
//     for (NodeID &v : clique) {
//       if (temp_status[v] == false) {
//         std::cout << "Overlap" << std::endl;
//         return;
//       }
//
//       temp_status[v] == false;
//       check_nodes--;
//     }
//
//     if (!safeClique(clique)) {
//       printVectorSet(clique);
//       std::cout << "Invalid clique" << std::endl;
//       return;
//     }
//   }
//
//   if (check_nodes > 0) {
//     std::cout << "Remaining nodes: " << check_nodes << std::endl;
//   }
// }
//



// void redu_structure::addKernelCliques(std::vector<std::vector<int>> &clique_set){
//
//   for (unsigned int i = 0; i < clique_set.size(); i++){
//       std::vector<NodeID> clique;
//
//       for (unsigned int j = 0; j < clique_set[i].size(); j++){
//           int v = clique_set[i][j];
//           NodeID old_v = new_to_old_map[v];
//
//           solve_node_clique[old_v] = false;
//           clique.push_back(old_v);
//       }
//       std::sort(clique.begin(), clique.end());
//
//       addCliqueToCover(clique);
//   }
// }
//
// void redu_structure::addCrownCliques(std::vector<std::vector<NodeID>> &crown_cliques, std::vector<std::vector<int>> &clique_set) {
//
//   for (unsigned int i = 0; i < clique_set.size(); i++){
//       std::vector<NodeID> clique;
//
//       for (unsigned int j = 0; j < clique_set[i].size(); j++){
//           int v = clique_set[i][j];
//           NodeID old_v = new_to_old_map[v];
//           // std::cout << old_v << std::endl;
//           // solve_node_clique[old_v] = false;
//           clique.push_back(old_v);
//       }
//       std::sort(clique.begin(), clique.end());
//
//       // printVectorSet(clique);
//       addClique(clique);
//       removeVertexSet(clique);
//
//       crown_cliques.push_back(clique);
//   }
//
// }
//
// unsigned int redu_structure::adj_size(NodeID v) {
//
//   unsigned int size = 0;
//   for (NodeID u: adj_list[v]) {
//     if (node_status[u]) { size++; }
//   }
//
//   return size;
// }
//
// std::vector<NodeID> redu_structure::curr_adj_list(NodeID v) {
//
//   std::vector<NodeID> curr_adj_list;
//
//   for (NodeID u : adj_list[v]){
//     if (node_status[u]) { curr_adj_list.push_back(u); }
//   }
//
//   return curr_adj_list;
// }
//
//
// void redu_structure::removeVertex(NodeID v) {
//   // removes vertex
//
//   node_status[v] = false;
//   remaining_nodes--;
//
//   if (!node_mis.empty() && node_mis[v]) curr_mis--;
// }
//
// void redu_structure::addVertex(NodeID v) {
//   // adds vertex
//
//   node_status[v] = true;
//   remaining_nodes++;
//
//   if (!node_mis.empty() && node_mis[v]) curr_mis++;
// }
//
// // void redu_structure::addClique(std::vector<NodeID> &clique) {
// //
// //   std::sort(clique.begin(), clique.end());
// //   clique_cover.push_back(clique);
// //   unsigned int cliqueID = clique_cover.size() - 1;
// //
// //   for (NodeID u : clique){
// //       node_clique[u] = cliqueID;
// //   }
// //
// // }
//
// void redu_structure::addClique(std::vector<NodeID> &clique) {
//
//   for (NodeID u : clique) {
//     node_clique[u] = next_cliqueID;
//   }
//
//   next_cliqueID++;
// }
//
// void redu_structure::addCliqueToCover(std::vector<NodeID> &clique) {
//   // adds clique so solve node structure
//
//   for (NodeID u : clique) {
//     solve_node_clique[u] = next_solvecliqueID;
//   }
//   clique_cover.push_back(clique);
//   next_solvecliqueID++;
//   // std::cout << next_solvecliqueID << std::endl;
// }
//
//
//
// // std::vector<NodeID> redu_structure::pop_clique() {
// //
// //   std::vector<NodeID> clique = clique_cover.back();
// //   clique_cover.pop_back();
// //   return clique;
// // }
//
// void redu_structure::pop_clique(std::vector<NodeID> &clique) {
//
//   for (NodeID u : clique) {
//     node_clique[u] = node_clique.size();
//   }
//
//   next_cliqueID--;
// }
//
// unsigned int redu_structure::getCliqueID(NodeID &v) {
//
//   return node_clique[v];
// }
//
// std::vector<NodeID> redu_structure::getClique(NodeID &v) {
//
//   unsigned int cliqueID = getCliqueID(v);
//   std::vector<NodeID> clique = clique_cover[cliqueID];
//
//   return clique;
// }
//
// void redu_structure::replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique){
//
//   std::sort(new_clique.begin(), new_clique.end());
//
//   for (NodeID a : new_clique) {
//     solve_node_clique[a] = cliqueID;
//   }
//   clique_cover[cliqueID] = new_clique;
// }
//
//
// void redu_structure::removeVertexSet(std::vector<NodeID> &S) {
//
//   for (NodeID v : S) { removeVertex(v); };
// }
//
// void redu_structure::addVertexSet(std::vector<NodeID> &S) {
//   for (NodeID v : S) { addVertex(v); };
// }
//
// // void redu_structure::clearScratch(std::vector<bool> &scratch) {
// //
// //   std::fill(scratch.begin(), scratch.end(), false);
// // }
//
void redu_structure::printAdjList() {

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
//
// void redu_structure::printAdjList(NodeID v) {
//
//   if (!node_status[v]) {
//     std::cout << v << " /notin adj_list" << std::endl;
//     return;
//   }
//
//   std::cout << "N(" << v << "): [";
//   for (NodeID u : adj_list[v]) {
//     if (!node_status[u]) { continue; }
//     std::cout << u << ", ";
//   }
//   std::cout << "]" << std::endl;
// }
//
// void redu_structure::printNeighborhood(NodeID v) {
//
//   printAdjList(v);
//   for (NodeID u : adj_list[v]) { std::cout << "  "; printAdjList(u); }
// }
//
// void redu_structure::printVectorSet(std::vector<NodeID> S){
//
//   for (unsigned int i = 0; i < S.size(); i++){
//     std::cout << S[i] << ", ";
//   }
//   std::cout << std::endl;
// }
