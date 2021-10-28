
#include <algorithm>
#include <iostream>
#include <fstream>

#include "dom_reduction.h"

bool dom_reduction::nodeDominates(redu_vcc &reduVCC, NodeID &v, NodeID &a){

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  unsigned int i = 0;
  unsigned int j = 0;

  while (j < adj_list[a].size()){
      if (!reduVCC.node_status[adj_list[a][j]] || adj_list[a][j] == v){
          j++;
          continue;
      }
      if (i == adj_list[v].size()){
          return false;
      }
      if (!reduVCC.node_status[adj_list[v][i]] || adj_list[v][i] == a){
          i++;
          continue;
      }
      if (adj_list[v][i] > adj_list[a][j]){
          return false;
      }
      if (adj_list[v][i] < adj_list[a][j]){
          i++;
      }
      else {
          i++;
          j++;
      }
  }

  return true;
}


bool dom_reduction::validDOM(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

    if (reduVCC.adj_size(v) <= 2) return false;

    for (NodeID a : adj_list[v]) {
      if (!reduVCC.node_status[a]) continue;
      if (reduVCC.adj_size(v) < reduVCC.adj_size(a)) continue;

      if (dom_reduction::nodeDominates(reduVCC, v, a)) {
        u = a;
        return true;
      }
    }
    return false;
}


void dom_reduction::reduce(  redu_vcc &reduVCC,
                           NodeID &node_v, NodeID &node_u ){

  type = "dom";

  v = node_v;
  u = node_u;

  deg = reduVCC.adj_size(v);

  reduVCC.removeVertex(v);
  reduVCC.fold_node[v] = true;

}

void dom_reduction::reduce(  redu_vcc &reduVCC, vertex_queue *queue,
                           NodeID &node_v, NodeID &node_u ){
  reduce( reduVCC, node_v, node_u);
  queue->adjust_queue(reduVCC, v);
}

void dom_reduction::unfold( redu_vcc &reduVCC){

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

////  assert(u < reduVCC.solve_node_clique.size());
  unsigned int fold_cliqueID = reduVCC.solve_node_clique[u];
  std::vector<NodeID> fold_clique = reduVCC.clique_cover[fold_cliqueID];

  fold_clique.push_back(v);
  reduVCC.replaceClique(fold_cliqueID, fold_clique);
}

void dom_reduction::unreduce( redu_vcc &reduVCC){
  reduVCC.addVertex(v);
  reduVCC.fold_node[v] = false;
}
