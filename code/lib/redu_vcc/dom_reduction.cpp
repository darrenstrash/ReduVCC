
#include <algorithm>
#include <iostream>
#include <fstream>

#include "dom_reduction.h"

bool dom_reduction::nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &a){

  unsigned int i = 0;
  unsigned int j = 0;

  while (j < adj_list[a].size()){
      if (adj_list[a][j] == v){
          j++;
          continue;
      }
      if (i == adj_list[v].size()){
          return false;
      }
      if (adj_list[v][i] == a){
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

    if (adj_list[v].size() <= 2) {return false;}

    for (NodeID a : adj_list[v]) {
      if (adj_list[v].size() < adj_list[a].size()) { continue; }
      if (dom_reduction::nodeDominates(adj_list, v, a)) {
        u = a;
        return true;
      }
    }

    return false;
}


void dom_reduction::reduce(graph_access &G,  redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;
  u = node_u;

  reduVCC.removeVertex(v);

}

void dom_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

  unsigned int cliqueID = reduVCC.getCliqueID(u);
  std::vector<NodeID> clique = reduVCC.getClique(u);

  clique.push_back(v);
  reduVCC.replaceClique(cliqueID, clique);
}
