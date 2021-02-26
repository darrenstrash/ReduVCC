
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

redu_vcc::redu_vcc(graph_access &G) {

  // produce adjacency list
  generateAdjList(G);
  // assign status of nodes
  node_status.assign(G.number_of_nodes(), true);
  fold_node.assign(G.number_of_nodes(), false);
  remaining_nodes = G.number_of_nodes();

  // allocate for graph cover
  node_clique.resize(G.number_of_nodes());

  // allocate two scratch vectors
  scratch1.assign(G.number_of_nodes(), false);
  scratch2.assign(G.number_of_nodes(), false);

  // assign first cliqueID to 0
  next_cliqueID = 0;
}

void redu_vcc::build_cover(graph_access &G){
  /* Constructs clique cover from node_clique mapping. */

  // clears clique_cover
  clique_cover.clear();
  clique_cover.assign(next_cliqueID, {});

  forall_nodes(G, v) {
    if (fold_node[v]) { continue; } // if node in fold, skip
    if (node_status[v]) { continue; } // if node still in graph, skip

    unsigned int cliqueID = node_clique[v];
    clique_cover[cliqueID].push_back(v);
  } endfor

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

void redu_vcc::assignMaps(graph_access &G) {

  old_to_new_map.clear();
  old_to_new_map.resize(G.number_of_nodes());
  new_to_old_map.clear();
  new_to_old_map.resize(G.number_of_nodes());

  int j = 0;
  for (unsigned int i = 0; i < G.number_of_nodes(); i++){
      if (!node_status[i]){
          continue;
      }
      old_to_new_map[i] = j;
      new_to_old_map[j] = i;
      j++;
  }
}

void redu_vcc::buildKernel(graph_access &G) {

  assignMaps(G);

  kernel_adj_list.clear();
  kernel_adj_list.resize(remaining_nodes);
  kernel_edges = 0;

  for (unsigned int i = 0; i < G.number_of_nodes(); i++) {
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
}

void redu_vcc::addVertex(NodeID v) {
  // adds vertex

  node_status[v] = true;
  remaining_nodes++;
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
