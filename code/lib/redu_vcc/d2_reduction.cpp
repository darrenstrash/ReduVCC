
#include <algorithm>
#include <iostream>
#include <fstream>

#include "d2_reduction.h"


void d2_reduction::assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w) {
  // assigns NodeID u to the neighbor of v with the LEAST neighbors

  NodeID x = adj_list[v][0];
  NodeID y = adj_list[v][1];

  if (adj_list[x].size() >= adj_list[y].size()) {u = y; w = x;}
  else {u = x; w = y;}
}

bool d2_reduction::isTriangle(std::vector<std::vector<NodeID>> &adj_list, NodeID &u, NodeID &w){

  for (NodeID x : adj_list[u]){
      if (x == w) {return true;}
  }
  return false;
}

bool d2_reduction::validNeighbors(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID & w){

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;
  std::vector<bool> &scratch2 = reduVCC.scratch2;

  // scratch1 = N[u] / N[w]
  // scratch 2 = N[w] / N[u]
  for (NodeID x : adj_list[u]) {scratch1[x] = true;}
  for (NodeID x : adj_list[w]) {scratch2[x] = true; scratch1[x] = false;}
  for (NodeID x : adj_list[u]) {scratch2[x] = false;}

  for (NodeID x : adj_list[u]){
      if (!scratch1[x]) {continue;}
      for (NodeID y : adj_list[x]) {
          if (scratch2[y]) {
            reduVCC.clearScratch(scratch1);
            reduVCC.clearScratch(scratch2);
            return false;
          }
      }
  }

  reduVCC.clearScratch(scratch1);
  reduVCC.clearScratch(scratch2);
  return true;
}


bool d2_reduction::validD2(redu_vcc &reduVCC, NodeID &v){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

    if (adj_list[v].size() != 2) {return false;}

    NodeID u; NodeID w;
    d2_reduction::assignNodes(adj_list, v, u, w);

    if (d2_reduction::isTriangle(adj_list, u, w)) {return false;}

    return reduction::uncrossedSets(reduVCC, u, w);

    //
    if (d2_reduction::validNeighbors(reduVCC, v, u, w)) {
    //     // std::cout << v << ": " << u << ", " << w << std::endl;
    //     for (NodeID x : adj_list[u]) {scratch1[x] = false;}
    //     for (NodeID x : adj_list[w]) {scratch2[x] = false;}
        return true;}
    //
    // for (NodeID x : adj_list[u]) {scratch1[x] = false;}
    // for (NodeID x : adj_list[w]) {scratch2[x] = false;}
    return false;

}

void d2_reduction::foldD2(redu_vcc &reduVCC) {

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
    std::vector<bool> &scratch1 = reduVCC.scratch1;

    for (NodeID x : adj_list[w]) {scratch1[x] = true;}

    for (NodeID x : adj_list[u]) {
        if (scratch1[x]) {continue;}
        adj_list[w].push_back(x);
        adj_list[x].push_back(w);
        std::sort(adj_list[x].begin(), adj_list[x].end());
    }
    std::sort(adj_list[w].begin(), adj_list[w].end());

    reduVCC.clearScratch(scratch1);
    reduVCC.removeVertex(u);

}

void d2_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;
  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  d2_reduction::assignNodes(adj_list, v, u, w);

  reduVCC.removeVertex(v);

  N_u = adj_list[u];
  foldD2(reduVCC);

}

void d2_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){
    // std::cout << "Unreducing D2... " << std::endl;

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
    std::vector<bool> &scratch1 = reduVCC.scratch1;

    unsigned int fold_cliqueID = reduVCC.getCliqueID(w);
    std::vector<NodeID> fold_clique = reduVCC.getClique(w);

    // std::cout << "N[u]: ";
    // reduVCC.printVectorSet(N_u);
    // std::cout << "fold clique: ";
    // reduVCC.printVectorSet(fold_clique);

    NodeID x = u;   // vertex x will be connected "externally" -- outside [v,u,w] or no connection
    NodeID y = w;   // vertex y will be connected to v

    // bool w_in = false;

    // tests to see if fold_clique \subset N_u,
    // if not then must be that fold_clique \subset N_w
    // if (!isSubset(reduVCC, fold_clique, N_u)) {
    //   std::cout << "switch" << std::endl;
    //   x = w; y = u;
    // }
    for (NodeID a : N_u) {scratch1[a] = true;}
    for (NodeID a : fold_clique) {
        if (a == w){
            continue;
        }
        if (!scratch1[a]) {
            // w_in = true;
            x = w; y = u;
        }
    }
    reduVCC.clearScratch(scratch1);


    for (unsigned int i = 0; i < fold_clique.size(); i++){
        if (fold_clique[i] == w){
            fold_clique[i] = x;
            reduVCC.node_clique[x] = fold_cliqueID;
            std::sort(fold_clique.begin(), fold_clique.end());
            break;
        }
    }
    reduVCC.replaceClique(fold_cliqueID, fold_clique);

    std::vector<NodeID> new_clique {v, y};
    // new_clique.push_back(v);
    // new_clique.push_back(y);
    reduVCC.addClique(new_clique);

    // std::cout << "new_clique: ";
    // reduVCC.printVectorSet(fold_clique);
    // std::cout << "new_clique2: ";
    // reduVCC.printVectorSet(new_clique);
    //
    // std::cout << std::endl;

}
