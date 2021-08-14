
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

  std::vector<NodeID> adj_2 = reduVCC.curr_adj_list(a);
std::vector<NodeID> adj_1  = reduVCC.curr_adj_list(v);
adj_2.push_back(a);
adj_1.push_back(v);
// std::sort(adj_list[p].begin(), adj_list[p].end());
// reduVCC.printVectorSet(adj_1);
// reduVCC.printVectorSet(adj_2);

  // if (!dom_reduction::isSubset(reduVCC, adj_2, adj_1)) std::cout << "error" << std::endl;
  //
  // reduVCC.printAdjList(v);
  // reduVCC.printAdjList(a);
  return true;
}


bool dom_reduction::validDOM(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

    // if (adj_list[v].size() <= 2) {return false;}
    if (reduVCC.adj_size(v) <= 2) return false;

    for (NodeID a : adj_list[v]) {
      if (!reduVCC.node_status[a]) continue;
      // if (adj_list[v].size() < adj_list[a].size()) { continue; }
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

  // std::cout << "success" << std::endl;
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
