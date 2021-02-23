
#include <algorithm>
#include <iostream>
#include <fstream>

#include "d2_reduction.h"


void d2_reduction::assignNodes(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID &w) {
  // assigns NodeID u to the neighbor of v with the LEAST neighbors

  std::vector<NodeID> N_v = reduVCC.curr_adj_list(v);
  NodeID x = N_v[0];
  NodeID y = N_v[1];

  if (reduVCC.adj_size(x) >= reduVCC.adj_size(y)) {u = y; w = x;}
  else {u = x; w = y;}
}

bool d2_reduction::isTriangle(redu_vcc &reduVCC, NodeID &u, NodeID &w){

  for (NodeID x : reduVCC.adj_list[u]){
      if (!reduVCC.node_status[x]) { continue; }
      if (x == w) {return true;}
  }
  return false;
}


bool d2_reduction::validD2(redu_vcc &reduVCC, NodeID &v){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

    if (reduVCC.adj_size(v) != 2) {return false;}

    NodeID u; NodeID w;
    d2_reduction::assignNodes(reduVCC, v, u, w);

    if (d2_reduction::isTriangle(reduVCC, u, w)) {return false;}

    return reduction::uncrossedSets(reduVCC, u, w);


}

void d2_reduction::foldD2(redu_vcc &reduVCC) {

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
    std::vector<bool> &scratch1 = reduVCC.scratch1;

    for (NodeID x : adj_list[w]) {scratch1[x] = true;}

    for (NodeID x : adj_list[u]) {
        if (!reduVCC.node_status[x]) { continue; }

        N_u.push_back(x);
        if (scratch1[x]) {
          continue;}
        disjoint.push_back(x);

        adj_list[w].push_back(x);
        adj_list[x].push_back(w);
        std::sort(adj_list[x].begin(), adj_list[x].end());
    }
    std::sort(adj_list[w].begin(), adj_list[w].end());

    reduVCC.clearScratch(scratch1);
    reduVCC.removeVertex(u);
    reduVCC.fold_node[u] = true;

}

void d2_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;
  num_folded_cliques++;
  // std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  // reduVCC.printNeighborhood(v);

  d2_reduction::assignNodes(reduVCC, v, u, w);
  // std::cout << u << ", " << w << std::endl;
  // std::cout << reduVCC.node_status[w] << std::endl;

  reduVCC.removeVertex(v);
  reduVCC.fold_node[v] = true;

  // N_u = reduVCC.curr_adj_list(u);
  foldD2(reduVCC);
  // reduVCC.printAdjList(w);
  // std::cout << std::endl;

}

void d2_reduction::unfold(graph_access &G, redu_vcc &reduVCC) {
  // std::cout << "Unreducing D2... " << std::endl;

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

  unsigned int fold_cliqueID = reduVCC.solve_node_clique[w];
  std::vector<NodeID> fold_clique = reduVCC.clique_cover[fold_cliqueID];

  NodeID x = u;   // vertex x will be connected "externally" -- outside [v,u,w] or no connection
  NodeID y = w;   // vertex y will be connected to v

  // bool w_in = false;

  // tests to see if fold_clique \subset N_u
  for (NodeID a : N_u) {scratch1[a] = true;}
  for (NodeID a : fold_clique) {
      if (a == w){
          continue;
      }
      // \exists a \in fold_clique, a \notin N_u
      if (!scratch1[a]) {
          // fold_clique \subset N_w
          x = w; y = u;
      }
  }
  reduVCC.clearScratch(scratch1);


  for (unsigned int i = 0; i < fold_clique.size(); i++){
      if (fold_clique[i] == w){
          fold_clique[i] = x;
          reduVCC.solve_node_clique[x] = fold_cliqueID;
          std::sort(fold_clique.begin(), fold_clique.end());
          break;
      }
  }
  // reduVCC.printVectorSet(fold_clique);
  reduVCC.replaceClique(fold_cliqueID, fold_clique);

  std::vector<NodeID> new_clique {v, y};
  reduVCC.addCliqueToCover(new_clique);
  // reduVCC.printVectorSet(new_clique);

}

void d2_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){

    for (NodeID a : disjoint) { reduVCC.scratch1[a] = true; }

    for (unsigned int i = 0; i < reduVCC.adj_list[w].size(); i++) {
      NodeID x = reduVCC.adj_list[w][i];
      if (reduVCC.scratch1[x]) {
        reduVCC.adj_list[w].erase(reduVCC.adj_list[w].begin() + i);
      }
      for (unsigned int j = 0; j < reduVCC.adj_list[x].size(); j++) {
        NodeID y = reduVCC.adj_list[x][j];
        if (y == x) {
          reduVCC.adj_list[x].erase(reduVCC.adj_list[x].begin() + j);
        }
      }
    }

    reduVCC.addVertex(v);
    reduVCC.fold_node[v] = false;
    reduVCC.addVertex(u);
    reduVCC.fold_node[u] = false;
}
