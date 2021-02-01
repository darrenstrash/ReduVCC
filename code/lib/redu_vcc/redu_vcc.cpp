
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

void redu_vcc::build(graph_access &G) {

  // produce adjacency list
  generateAdjList(G);
  // assign status of nodes
  node_status.assign(G.number_of_nodes(), true);
  remaining_nodes = G.number_of_nodes();

  // allocate for graph cover
  node_clique.resize(G.number_of_nodes());

  scratch1.assign(G.number_of_nodes(), false);
  scratch2.assign(G.number_of_nodes(), false);

  // printAdjList();
}

bool redu_vcc::cliqueInG(graph_access &G, std::vector<NodeID> &clique) {

  for (unsigned int i = 0; i < clique.size() -1; i++) {
    NodeID v = clique[i];
    unsigned int j = i + 1;
    forall_out_edges(G, e, v){
        NodeID u = G.getEdgeTarget(e);

        // std::cout << u << ", " << clique[j] << std::endl;
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

  std::cout << clique_cover.size() << std::endl;
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

    if (!cliqueInG(G, clique)) { std::cout << "Invalid clique" << std::endl; return; }
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

          removeVertex(old_v);
          clique.push_back(old_v);
      }
      std::sort(clique.begin(), clique.end());

      addClique(clique);
  }
}


void redu_vcc::removeVertex(NodeID v) {
  // removes vertex from adjacency list

  node_status[v] = false;
  for (NodeID u : adj_list[v]){
      for (unsigned int i = 0; i < adj_list[u].size(); i++){
          NodeID z = adj_list[u][i];
          if (z == v){
              adj_list[u].erase(adj_list[u].begin() + i);
              continue;
          }
      }
  }
  adj_list[v].erase(adj_list[v].begin(), adj_list[v].end());

  remaining_nodes--;
}

void redu_vcc::addClique(std::vector<NodeID> &clique) {

  std::sort(clique.begin(), clique.end());
  clique_cover.push_back(clique);
  unsigned int cliqueID = clique_cover.size() - 1;

  for (NodeID u : clique){
      node_clique[u] = cliqueID;
  }

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
    node_clique[a] = cliqueID;
  }
  clique_cover[cliqueID] = new_clique;
}


void redu_vcc::removeVertexSet(std::vector<NodeID> &S) {

  for (NodeID v : S) { removeVertex(v); };
}

void redu_vcc::clearScratch(std::vector<bool> &scratch) {

  std::fill(scratch.begin(), scratch.end(), false);
}

void redu_vcc::printAdjList() {

  for (unsigned int i = 0; i < node_status.size(); i++) {
    if (!node_status[i]) { continue; }
    std::cout << "N(" << i << "): [";
    for (NodeID u : adj_list[i]) {
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
