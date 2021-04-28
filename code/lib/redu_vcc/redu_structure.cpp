
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

// void redu_structure::generateAdjList(graph_access &G) {
//   /* Generates adjacency list from graph */
//
//   forall_nodes(G, v){
//       std::vector<NodeID> N_v;    // open neighborhood of v
//
//       forall_out_edges(G, e, v){
//           NodeID u = G.getEdgeTarget(e);
//           N_v.push_back(u);
//
//       } endfor
//       adj_list.push_back(N_v);
//
//   } endfor
// }

// redu_structure::redu_structure(graph_access &G) {
//
//   // produce adjacency list
//   generateAdjList(G);
//   // assign status of nodes
//   num_nodes = G.number_of_nodes();
//   node_status.assign(num_nodes, true);
//   fold_node.assign(num_nodes, false);
//   remaining_nodes = num_nodes;
//   // allocate for graph cover
//   node_clique.resize(num_nodes);
//
//   // initialize mis mapping to 0
//   curr_mis = 0;
//
//   // allocate two scratch vectors
//   scratch1.assign(num_nodes, false);
//   scratch2.assign(num_nodes, false);
//
//   // assign first cliqueID to 0
//   next_cliqueID = 0;
// }

void redu_structure::init() {

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

void redu_structure::build_cover(graph_access &G){
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

bool redu_structure::cliqueInG(graph_access &G, std::vector<NodeID> &clique) {
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

void redu_structure::validateCover(graph_access &G) {

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


unsigned int redu_structure::adj_size(NodeID v) {

  unsigned int size = 0;
  for (NodeID u: adj_list[v]) {
    if (node_status[u]) { size++; }
  }

  return size;
}

std::vector<NodeID> redu_structure::curr_adj_list(NodeID v) {

  std::vector<NodeID> curr_adj_list;

  for (NodeID u : adj_list[v]){
    if (node_status[u]) { curr_adj_list.push_back(u); }
  }

  return curr_adj_list;
}


void redu_structure::removeVertex(NodeID v) {
  // removes vertex

  node_status[v] = false;
  remaining_nodes--;

  if (!node_mis.empty() && node_mis[v]) curr_mis--;
}

void redu_structure::addVertex(NodeID v) {
  // adds vertex

  node_status[v] = true;
  remaining_nodes++;

  if (!node_mis.empty() && node_mis[v]) curr_mis++;
}

// void redu_structure::addClique(std::vector<NodeID> &clique) {
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

void redu_structure::addClique(std::vector<NodeID> &clique) {

  for (NodeID u : clique) {
    node_clique[u] = next_cliqueID;
  }

  next_cliqueID++;
}

void redu_structure::addCliqueToCover(std::vector<NodeID> &clique) {
  // adds clique so solve node structure

  for (NodeID u : clique) {
    solve_node_clique[u] = next_solvecliqueID;
  }
  clique_cover.push_back(clique);
  next_solvecliqueID++;
  // std::cout << next_solvecliqueID << std::endl;
}



// std::vector<NodeID> redu_structure::pop_clique() {
//
//   std::vector<NodeID> clique = clique_cover.back();
//   clique_cover.pop_back();
//   return clique;
// }

void redu_structure::pop_clique(std::vector<NodeID> &clique) {

  for (NodeID u : clique) {
    node_clique[u] = node_clique.size();
  }

  next_cliqueID--;
}

unsigned int redu_structure::getCliqueID(NodeID &v) {

  return node_clique[v];
}

std::vector<NodeID> redu_structure::getClique(NodeID &v) {

  unsigned int cliqueID = getCliqueID(v);
  std::vector<NodeID> clique = clique_cover[cliqueID];

  return clique;
}

void redu_structure::replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique){

  std::sort(new_clique.begin(), new_clique.end());

  for (NodeID a : new_clique) {
    solve_node_clique[a] = cliqueID;
  }
  clique_cover[cliqueID] = new_clique;
}


void redu_structure::removeVertexSet(std::vector<NodeID> &S) {

  for (NodeID v : S) { removeVertex(v); };
}

void redu_structure::addVertexSet(std::vector<NodeID> &S) {
  for (NodeID v : S) { addVertex(v); };
}

// void redu_structure::clearScratch(std::vector<bool> &scratch) {
//
//   std::fill(scratch.begin(), scratch.end(), false);
// }

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

void redu_structure::printAdjList(NodeID v) {

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

void redu_structure::printNeighborhood(NodeID v) {

  printAdjList(v);
  for (NodeID u : adj_list[v]) { std::cout << "  "; printAdjList(u); }
}

void redu_structure::printVectorSet(std::vector<NodeID> S){

  for (unsigned int i = 0; i < S.size(); i++){
    std::cout << S[i] << ", ";
  }
  std::cout << std::endl;
}
