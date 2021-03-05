
#include <algorithm>
#include <iostream>
#include <fstream>

#include "twin_reduction.h"


void twin_reduction::assignNodes(redu_vcc &reduVCC, NodeID &v, NodeID &w, NodeID &x, NodeID &y) {
  // assigns NodeID w to the neighbor of v with the LEAST neighbors

std::vector<NodeID> N_v = reduVCC.curr_adj_list(v);

  NodeID a = N_v[0];
  NodeID b = N_v[1];
  NodeID c = N_v[2];

  std::vector<NodeID> neighbors {a, b, c};
  std::vector<unsigned int> neighbor_sizes {reduVCC.adj_size(a), reduVCC.adj_size(b), reduVCC.adj_size(c)};
  std::vector<unsigned int> neighbor_indices { 0, 1, 2 };

  std::sort(neighbor_indices.begin(), neighbor_indices.end(),
    [neighbor_sizes](unsigned int i, unsigned int j) {

        unsigned int s1 = neighbor_sizes[i];
        unsigned int s2 = neighbor_sizes[j];
        return s1 < s2;
    });

    w = neighbors[neighbor_indices[0]];
    x = neighbors[neighbor_indices[1]];
    y = neighbors[neighbor_indices[2]];

}


bool twin_reduction::validNeighbors(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y){

  for (NodeID a : reduVCC.adj_list[w]) {
    if (!reduVCC.node_status[a]) continue;
    if (a == y || a == x) return true;
  }
  for (NodeID a : reduVCC.adj_list[x]) {
    if (!reduVCC.node_status[a]) continue;
    if (a == y) { return true; }
  }

  return reduction::uncrossedSets(reduVCC, w, x) && reduction::uncrossedSets(reduVCC, w, y) && reduction::uncrossedSets(reduVCC, x, y);
}


bool twin_reduction::twinFound( redu_vcc &reduVCC,
                                NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y) {

      std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

      for (NodeID p : adj_list[w]){
      if (!reduVCC.node_status[p]) continue;
      if (p == v) continue;

      // if (adj_list[p].size() != 3) continue;
      if (reduVCC.adj_size(p) != 3) continue;

      for (NodeID q : adj_list[x]){
          if (!reduVCC.node_status[q]) continue;
          if (q > p) break;

          if (q == p){
              for (NodeID r : adj_list[y]){
                  if (!reduVCC.node_status[r]) continue;
                  if (r > p) break;

                  if (r == p){
                      u = p;
                      // std::cout << "fold found" << std::endl;
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

    // if (adj_list[v].size() != 3) return false;
    if (reduVCC.adj_size(v) != 3) return false;

    NodeID w;
    NodeID x;
    NodeID y;
    twin_reduction::assignNodes(reduVCC, v, w, x, y);

    // return false;

    return twin_reduction::twinFound(reduVCC, v, u, w, x, y) && twin_reduction::validNeighbors(reduVCC, v, u, w, x, y);
}


bool twin_reduction::removeType (redu_vcc &reduVCC) {

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  for (NodeID a : adj_list[w]) {
    if (!reduVCC.node_status[a]) continue;
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
    if (!reduVCC.node_status[a]) continue;
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
  reduVCC.fold_node[v] = true;
  reduVCC.removeVertex(u);
  reduVCC.fold_node[u] = true;

  // N_w = adj_list[w];
  // N_x = adj_list[x];

  merge_neighborhoods(reduVCC, disjoint, N_x, y, x);
  merge_neighborhoods(reduVCC, disjoint, N_w, y, w);

  reduVCC.removeVertex(w);
  reduVCC.fold_node[w] = true;
  reduVCC.removeVertex(x);
  reduVCC.fold_node[x] = true;

}

void twin_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;
  u = node_u;

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  twin_reduction::assignNodes(reduVCC, v, w, x, y);

  if (removeType(reduVCC)) {
    remove_type = true;
    removeTWIN(reduVCC);
    num_cliques += 2;
    return;
  }

  remove_type = false;
  foldTWIN(reduVCC);
  num_folded_cliques +=2;

}

void twin_reduction::unfoldTWIN(redu_vcc &reduVCC,
                                std::vector<NodeID> &partial_clique,
                                unsigned &cliqueID,
                                NodeID &a, NodeID &b, NodeID &c) {

  partial_clique.push_back(a);
  reduVCC.replaceClique(cliqueID, partial_clique);

  std::vector<NodeID> new_clique1 {v, b};
  reduVCC.addCliqueToCover(new_clique1);
  std::vector<NodeID> new_clique2 {u, c};
  reduVCC.addCliqueToCover(new_clique2);
}

void twin_reduction::unfold(graph_access &G, redu_vcc &reduVCC){
  if (remove_type) { return; }

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

  // unsigned int fold_cliqueID = reduVCC.getCliqueID(y);
  unsigned int fold_cliqueID = reduVCC.solve_node_clique[y];
  // std::vector<NodeID> fold_clique = reduVCC.getClique(y);
  std::vector<NodeID> fold_clique = reduVCC.clique_cover[fold_cliqueID];

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

void twin_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){

    reduVCC.addVertex(v);
    reduVCC.fold_node[v] = false;
    reduVCC.addVertex(u);
    reduVCC.fold_node[u] = false;
    reduVCC.addVertex(w);
    reduVCC.fold_node[w] = false;
    reduVCC.addVertex(x);
    reduVCC.fold_node[x] = false;

    if (remove_type) {
      reduVCC.addVertex(y);
      reduVCC.fold_node[y] = false;
      return;
    }


    for (NodeID a : disjoint) {
      for (unsigned int i = 0; i < reduVCC.adj_list[y].size(); i++) {
        if (reduVCC.adj_list[y][i] == a) {
          reduVCC.adj_list[y].erase(reduVCC.adj_list[y].begin() + i);
        }
      }
      for (unsigned int i = 0; i < reduVCC.adj_list[a].size(); i++) {
        if (reduVCC.adj_list[a][i] == y) {
          reduVCC.adj_list[a].erase(reduVCC.adj_list[a].begin() + i);
        }
      }
    }


}
