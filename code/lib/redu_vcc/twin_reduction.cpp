
#include <algorithm>
#include <iostream>
#include <fstream>

#include "twin_reduction.h"


void twin_reduction::assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &w, NodeID &x, NodeID &y) {
  // assigns NodeID u to the neighbor of v with the LEAST neighbors

  NodeID a = adj_list[v][0];
  NodeID b = adj_list[v][1];
  NodeID c = adj_list[v][2];

  std::vector<NodeID> neighbors {a, b, c};

  std::sort(neighbors.begin(), neighbors.end(),
    [adj_list](NodeID n1, NodeID n2) {
        std::vector<NodeID> adj1 = adj_list[n1];
        std::vector<NodeID> adj2 = adj_list[n2];
        return adj1.size() < adj2.size();
    });

    w = neighbors[0];
    x = neighbors[1];
    y = neighbors[2];

}


bool twin_reduction::validNeighbors(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y){

  for (NodeID a : reduVCC.adj_list[w]) {
    if (a == y || a == x) { std::cout << "remove found" << std::endl; return true; }
  }
  for (NodeID a : reduVCC.adj_list[x]) {
    if (a == y) { std::cout << "remove found" << std::endl; return true; }
  }

  return reduction::uncrossedSets(reduVCC, w, x) && reduction::uncrossedSets(reduVCC, w, y) && reduction::uncrossedSets(reduVCC, x, y);
}


bool twin_reduction::twinFound( std::vector<std::vector<NodeID>> &adj_list,
                                NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y) {

      for (NodeID p : adj_list[w]){
      if (p == v){
          continue;
      }
      if (adj_list[p].size() != 3){
          continue;
      }

      for (NodeID q : adj_list[x]){
          if (q > p){
              break;
          }
          if (q == p){
              for (NodeID r : adj_list[y]){
                  if (r > p){
                      break;
                  }
                  if (r == p){
                      u = p;
                      std::cout << "fold found" << std::endl;
                      return true;
                  }
              }
              break;
          }
      }
    }
    return false;
}


bool twin_reduction::validTWIN(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

    if (adj_list[v].size() != 3) {return false;}

    NodeID w;
    NodeID x;
    NodeID y;
    twin_reduction::assignNodes(adj_list, v, w, x, y);

    return twin_reduction::twinFound(adj_list, v, u, w, x, y) && twin_reduction::validNeighbors(reduVCC, v, u, w, x, y);
}


bool twin_reduction::removeType (std::vector<std::vector<NodeID>> &adj_list) {

  for (NodeID a : adj_list[w]) {
    if ( a == x) {
      edge_nodes = {w, x};
      nonedge_node = y;
      return true;
    }
    if ( a == y ) {
      edge_nodes = {w, y};
      nonedge_node = x;
      return true;
    }
  }

  for (NodeID a : adj_list[x]) {
    if ( a == y ) {
      edge_nodes = {x, y};
      nonedge_node = w;
      return true;
    }
  }

  return false;
}

void twin_reduction::removeTWIN (redu_vcc &reduVCC) {

  std::vector<NodeID> clique1 {v, edge_nodes[0], edge_nodes[1]};
  std::vector<NodeID> clique2 {u, nonedge_node};

  reduVCC.addClique(clique1);
  reduVCC.addClique(clique2);

  reduVCC.removeVertexSet(clique1);
  reduVCC.removeVertexSet(clique2);
}

void twin_reduction::foldTWIN(redu_vcc &reduVCC) {

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  reduVCC.removeVertex(v);
  reduVCC.removeVertex(u);

  N_w = adj_list[w];
  N_x = adj_list[x];

  merge_neighborhoods(reduVCC, y, x);
  merge_neighborhoods(reduVCC, y, w);

  reduVCC.removeVertex(w);
  reduVCC.removeVertex(y);

}

void twin_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;
  u = node_u;

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  twin_reduction::assignNodes(adj_list, v, w, x, y);

  if (removeType(adj_list)) {
    remove_type = true;
    removeTWIN(reduVCC);
    return;
  }

  remove_type = false;
  foldTWIN(reduVCC);

}

void twin_reduction::unfoldTWIN(redu_vcc &reduVCC,
                                std::vector<NodeID> &partial_clique,
                                unsigned &cliqueID,
                                NodeID &a, NodeID &b, NodeID &c) {

  partial_clique.push_back(a);
  reduVCC.replaceClique(cliqueID, partial_clique);

  std::vector<NodeID> new_clique1 {v, b};
  reduVCC.addClique(new_clique1);
  std::vector<NodeID> new_clique2 {u, c};
  reduVCC.addClique(new_clique2);
}

void twin_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){
  if (remove_type) { return; }

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

  unsigned int fold_cliqueID = reduVCC.getCliqueID(y);
  std::vector<NodeID> fold_clique = reduVCC.getClique(y);

  std::vector<NodeID> partial_clique; // fold_clique \setminus y
  for (NodeID a : fold_clique) {
    if ( a == y ) { continue; }
    partial_clique.push_back(a);
  }

  if (isSubset(reduVCC, partial_clique, N_w)) {
    unfoldTWIN(reduVCC, partial_clique, fold_cliqueID, w, x, y);
  }
  else if (isSubset(reduVCC, partial_clique, N_x)) {
    unfoldTWIN(reduVCC, partial_clique, fold_cliqueID, x, w, y);
  }
  else {
    unfoldTWIN(reduVCC, partial_clique, fold_cliqueID, y, w, x);
  }

}
