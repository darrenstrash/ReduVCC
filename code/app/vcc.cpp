/*****A*************************************************************************
 * vcc.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

#include "mis/kernel/branch_and_reduce_algorithm.h"

#include "ccp/Chalupa/cli.h"
#include "ccp/Chalupa/random_generator.h"
#include <time.h>


class Reduction{
private:
  std::vector<bool> scratch;
    
  std::vector<int> reductions;
    
  std::vector<NodeID> d2_fold;
  std::vector<std::vector<NodeID>> d2_folded_nodes;
  std::vector<std::vector<NodeID>> d2_folded_neighborhood;
  
  std::vector<NodeID> twin_fold;
  std::vector<std::vector<NodeID>> twin_folded_nodes;
  std::vector<std::vector<std::vector<NodeID>>> twin_folded_neighborhoods;
  
  std::vector<NodeID> dominated_node;
  std::vector<NodeID> dominating_node;
  
  std::vector<int> old_new_map;
  std::vector<NodeID> new_old_map;
  std::vector<std::vector<int>> int_adj_list;
  
  void removeVertex(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
  void addClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<NodeID> &clique);
  
  void reduceIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  int getIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  bool checkNeighbor(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u);
  void addIsolatedClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
  void removeIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
  
  void reduceDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  int getDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  bool checkValidDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
  void foldDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
  
  void reduceTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  int getTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  bool checkValidTwinFold(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c);
  void organizeTwinElements(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &w, NodeID &y, NodeID &z, NodeID &a, NodeID &b, NodeID &c);
  bool foundTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c);
  bool removeTypeTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &a, NodeID &b, NodeID &c);
  void removeTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c);
  void foldTwinVertices(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c);
  
  void reduceDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  int getDominant(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  bool checkDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u);
  void removeDominant(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u);
  
  void reduceCrown(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  void addCrownCliques(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<int>> &crown_cliques);
  void makeNewAdj(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  
  bool checkCliqueStatus(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<NodeID> &clique, NodeID &x, std::vector<NodeID> &n_v);
  
  void unreduceDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  void adjustDegreeTwoClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &x, std::vector<NodeID> &fold, std::vector<NodeID> &n_u);
  
  void unreduceTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
  void adjustTwinClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &c, std::vector<NodeID> &fold, std::vector<std::vector<NodeID>> &neighborhoods);
    
  void unreduceDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
    
public:
  std::vector<bool> node_status;
  std::vector<std::vector<NodeID>> clique_cover;
  std::vector<unsigned int> node_clique;
  
  Reduction(graph_access &G);
  
  int getGraphSize(graph_access &G);
  void makeSubGraph(graph_access&G, std::vector<std::vector<NodeID>> &adj_list, PartitionConfig &partition_config);
  void reduceGraph(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, PartitionConfig &partition_config);
  void unreduceGraph(graph_access &G, std::vector<std::vector<NodeID>> &adj_list);
};


Reduction::Reduction(graph_access &G){
  node_status.assign(G.number_of_nodes(), true);
  scratch.assign(G.number_of_nodes(), false);
  node_clique.resize(G.number_of_nodes());
  
  //adj_list_map.resize(G.number_of_nodes());
}


bool Reduction::checkCliqueStatus(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<NodeID> &clique, NodeID &x, std::vector<NodeID> &n_v){
    
  unsigned int i = 0;
  unsigned int j = 0;
  
  while (i < clique.size()){
    if (clique[i] == x){
      i++;
      continue;
    }
    if (j == n_v.size()){
      return false;
    }
    if (clique[i] < n_v[j]){
      return false;
    }
    if (clique[i] == n_v[j]){
      i++;
      j++;
    }
    else{
      j++;
    }
  }
  return true;
}


void Reduction::unreduceDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
  NodeID v = dominating_node.back();
  dominating_node.pop_back();
  
  NodeID u = dominated_node.back();
  dominated_node.pop_back();
  
  unsigned int cliqueID = node_clique[u];
  clique_cover[cliqueID].push_back(v);
  node_clique[v] = cliqueID;
  sort(clique_cover[cliqueID].begin(), clique_cover[cliqueID].end());
}


void Reduction::adjustTwinClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &c, std::vector<NodeID> &fold, std::vector<std::vector<NodeID>> &neighborhoods){
  
  unsigned int cliqueID = node_clique[c];
  std::vector<NodeID> clique = clique_cover[cliqueID];
  
  NodeID v = fold[0];
  NodeID u = fold[1];
  NodeID a = fold[2];
  NodeID b = fold[3];
  
  NodeID x;
  NodeID y;
  NodeID z;
  
  std::vector<NodeID> n_a = neighborhoods[0];
  std::vector<NodeID> n_b = neighborhoods[1];
    
  if (checkCliqueStatus(G, adj_list, clique, c, n_a)){
    x = a;
    y = b;
    z = c;
  }
  else if (checkCliqueStatus(G, adj_list, clique, c, n_b)){
    x = b;
    y = a;
    z = c;
  }
  else{
    x = c;
    y = a;
    z = b;
  }
  
  for (unsigned int i = 0; i < clique.size(); i++){
    if (clique[i] == c){
      clique[i] = x;
      node_clique[x] = cliqueID;
      std::sort(clique.begin(), clique.end());
      break;
    }
  }
  
  clique_cover[cliqueID] = clique;
    
  std::vector<NodeID> new_clique1;
  new_clique1.push_back(v);
  new_clique1.push_back(y);
  std::sort(new_clique1.begin(), new_clique1.end());
  clique_cover.push_back(new_clique1);
  node_clique[v] = clique_cover.size() - 1;
  node_clique[y] = clique_cover.size() - 1;
  
  std::vector<NodeID> new_clique2;
  new_clique2.push_back(u);
  new_clique2.push_back(z);
  std::sort(new_clique2.begin(), new_clique2.end());
  clique_cover.push_back(new_clique2);
  node_clique[u] = clique_cover.size() - 1;
  node_clique[z] = clique_cover.size() - 1;
  
}


void Reduction::unreduceTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    NodeID c = twin_fold.back();
    twin_fold.pop_back();
    
    std::vector<NodeID> fold = twin_folded_nodes.back();
    twin_folded_nodes.pop_back();
    
    std::vector<std::vector<NodeID>> neighborhoods = twin_folded_neighborhoods.back();
    twin_folded_neighborhoods.pop_back();
    
    adjustTwinClique(G, adj_list, c, fold, neighborhoods);
}


void Reduction::adjustDegreeTwoClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &x, std::vector<NodeID> &fold, std::vector<NodeID> &n_u){
    
    unsigned int cliqueID = node_clique[x];
    std::vector<NodeID> clique = clique_cover[cliqueID];
    
    NodeID u = fold[0];
    NodeID v = fold[1];
    NodeID w = fold[2];
    
    NodeID z = u;
    NodeID y = w;
    
//    unsigned int i = 0;
//    unsigned int j = 0;
//
//    while (i < clique.size()){
//        if (clique[i] == x){
//            i++;
//            continue;
//        }
//        if (j == n_u.size()){
//            z = w;
//            y = u;
//            break;
//        }
//        if (clique[i] < n_u[j]){
//            z = w;
//            y = u;
//            break;
//        }
//        if (clique[i] == n_u[j]){
//            i++;
//            j++;
//        }
//        else{
//            j++;
//        }
//    }
    
    if (!checkCliqueStatus(G, adj_list, clique, x, n_u)){
        z = w;
        y = u;
    }
    
    for (unsigned int i = 0; i < clique.size(); i++){
        if (clique[i] == x){
            clique[i] = z;
            node_clique[z] = cliqueID;
            std::sort(clique.begin(), clique.end());
            break;
        }
    }
    
    clique_cover[cliqueID] = clique;
    
    std::vector<NodeID> new_clique;
    new_clique.push_back(v);
    new_clique.push_back(y);
    std::sort(new_clique.begin(), new_clique.end());
    clique_cover.push_back(new_clique);
    node_clique[v] = clique_cover.size() - 1;
    node_clique[y] = clique_cover.size() - 1;
    
}


void Reduction::unreduceDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    NodeID x = d2_fold.back();
    d2_fold.pop_back();
    
    std::vector<NodeID> fold = d2_folded_nodes.back();
    d2_folded_nodes.pop_back();
    
    std::vector<NodeID> n_u = d2_folded_neighborhood.back();
    d2_folded_neighborhood.pop_back();
    
    adjustDegreeTwoClique(G, adj_list, x, fold, n_u);
}


void Reduction::unreduceGraph(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    for (int i = reductions.size() - 1; i >= 0; i--){
        if (reductions[i] == 1){
            unreduceDegreeTwo(G, adj_list);
        }
        else if (reductions[i] == 2){
            unreduceTwin(G, adj_list);
        }
        else if (reductions[i] == 3)
            unreduceDominance(G, adj_list);
    }
}


void Reduction::removeVertex(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v){
    
    node_status[v] = false;
    for (unsigned int i = 0; i < adj_list[v].size(); i++){
        NodeID u = adj_list[v][i];
        for (unsigned int j = 0; j < adj_list[u].size(); j++){
            NodeID z = adj_list[u][j];
            if (z == v){
                adj_list[u].erase(adj_list[u].begin() + j);
                continue;
            }
        }
    }
    adj_list[v].erase(adj_list[v].begin(), adj_list[v].end());
}


void Reduction::addClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<NodeID> &clique){
    
    clique_cover.push_back(clique);
    unsigned int cliqueID = clique_cover.size() - 1;
    for (unsigned int i = 0; i < clique.size(); i++){
        NodeID v = clique[i];
//        std::cout << v << " ";
        node_clique[v] = cliqueID;
        removeVertex(G, adj_list, v);
    }
//    std::cout << std::endl;
}


void Reduction::makeNewAdj(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
  old_new_map.clear();
  old_new_map.resize(G.number_of_nodes());
  new_old_map.clear();
  new_old_map.resize(G.number_of_nodes());

  int_adj_list.clear();
  
  int u = 0;
    
  for (unsigned int i = 0; i < G.number_of_nodes(); i++){
    if (!node_status[i]){
      continue;
    }
    old_new_map[i] = u;
    new_old_map[u] = i;
    u++;
  }
    
  int_adj_list.resize(u);
    
  for (unsigned int i = 0; i < G.number_of_nodes(); i++){
    if (!node_status[i]){
      continue;
    }
    int new_v = old_new_map[i];
    std::vector<int> adj;
    for (unsigned int j = 0; j < adj_list[i].size(); j++){
      NodeID u  = adj_list[i][j];
      int new_u = old_new_map[u];
      adj.push_back(new_u);
    }
    std::sort(adj.begin(), adj.end());
    int_adj_list[new_v] = adj;
  }
}


void Reduction::addCrownCliques(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<int>> &crown_cliques){

  for (unsigned int i = 0; i < crown_cliques.size(); i++){
    std::vector<NodeID> clique;

    for (unsigned int j = 0; j < crown_cliques[i].size(); j++){
      int v = crown_cliques[i][j];
      NodeID old_v = new_old_map[v];

      clique.push_back(old_v);
    }
    std::sort(clique.begin(), clique.end());

    addClique(G, adj_list, clique);
  }
}
  

void Reduction::reduceCrown(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){

  makeNewAdj(G, adj_list);

  // std::cout << int_adj_list.size() << std::endl;
  branch_and_reduce_algorithm b_and_r(int_adj_list, int_adj_list.size());
  if (b_and_r.lpCrownReduction()){
    addCrownCliques(G, adj_list, b_and_r.crown_cliques);
  }
  // std::cout << something << std::endl;
  // std::cout << b_and_r.rn << std::endl;
  // std::cout << " size " << int_adj_list.size() << std::endl;
  // // for (unsigned int i = 0; i < b_and_r.crown_cliques.size(); i++){
  // //   std::cout << i << " ";
  // // }
  // // std::cout << std::endl;
  // std::cout << G.number_of_nodes() - b_and_r.rn << std::endl;

  // std::cout << int_adj_list.size() << std::endl;
  
}


void Reduction::removeDominant(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u) {
    
    dominating_node.push_back(v);
    dominated_node.push_back(u);
    
    removeVertex(G, adj_list, v);
}


bool Reduction::checkDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u){
    
    unsigned int i = 0;
    unsigned int j = 0;
    
    while (j < adj_list[u].size()){
        if (adj_list[u][j] == v){
            j++;
            continue;
        }
        if (i == adj_list[v].size()){
            return false;
        }
        if (adj_list[v][i] == u){
            i++;
            continue;
        }
        if (adj_list[v][i] > adj_list[u][j]){
            return false;
        }
        if (adj_list[v][i] < adj_list[u][j]){
            i++;
        }
        else {
            i++;
            j++;
        }
    }
    return true;
}


int Reduction::getDominant(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    int num_dominant = 0;
    forall_nodes(G, v) {
        if (!node_status[v]){
            continue;
        }
        if (adj_list[v].size() >= 3){
            for (unsigned int i = 0; i < adj_list[v].size(); i++){
                NodeID u = adj_list[v][i];
                if (adj_list[v].size() < adj_list[u].size()){
                    continue;
                }
                if (checkDominance(G, adj_list, v, u)){
//                    std::cout << "found: " << v << ", " << u << std::endl;
//                    for (unsigned int j = 0; j < adj_list[v].size(); j++){
//                        std::cout << adj_list[v][j] << ", ";
//                    }
//                    std::cout << std::endl;
//                    for (unsigned int j = 0; j < adj_list[u].size(); j++){
//                        std::cout << adj_list[u][j] << ", ";
//                    }
//                    std::cout << std::endl;
                    removeDominant(G, adj_list, v, u);
                    num_dominant++;
                    reductions.push_back(3);
                    break;
                }
            }
        }
    } endfor
    return num_dominant;
}


void Reduction::reduceDominance(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    bool dominance = true;

    while (dominance){
        int num_dominant = getDominant(G, adj_list);
        if (num_dominant == 0){
            dominance = false;
        }
    }
//        getDominant(G, adj_list);
}


void Reduction::foldTwinVertices(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c){
    
    // NodeID w = adj_list[v][0];
    // NodeID y = adj_list[v][1];
    // NodeID z = adj_list[v][2];
    
    std::vector<NodeID> fold;
    fold.push_back(v);
    fold.push_back(u);
    fold.push_back(a);
    fold.push_back(b);
    fold.push_back(c);
    twin_folded_nodes.push_back(fold);
    //
    std::vector<std::vector<NodeID>> neighborhoods;
    neighborhoods.push_back(adj_list[a]);
    neighborhoods.push_back(adj_list[b]);
    //    neighborhoods.push_back(adj_list[c]);
    twin_folded_neighborhoods.push_back(neighborhoods);
    
    removeVertex(G, adj_list, v);
    removeVertex(G, adj_list, u);
    
//    std::cout << "here" << std::endl;
    for (unsigned int i = 0; i < adj_list[c].size(); i++){
        NodeID r  = adj_list[c][i];
        scratch[r] = true;
    }
//    std::cout << "here2" << std::endl;
    for (unsigned int i = 0; i < adj_list[b].size(); i++){
        NodeID q = adj_list[b][i];
        
        if (q == c){
            continue;
        }
        if (scratch[q]){
            continue;
        }
        adj_list[c].push_back(q);
        adj_list[q].push_back(c);
        
        std::sort(adj_list[q].begin(), adj_list[q].end());
    }
//    std::cout << "here3" << std::endl;
    for (unsigned int i = 0; i < adj_list[c].size(); i++){
        NodeID r  = adj_list[c][i];
        scratch[r] = true;
    }
//    std::cout << "here4" << std::endl;
    for (unsigned int i = 0; i < adj_list[a].size(); i++){
        NodeID p = adj_list[a][i];
        
        if (p == c){
            continue;
        }
        if (scratch[p]){
            continue;
        }
        adj_list[c].push_back(p);
        adj_list[p].push_back(c);
        
        std::sort(adj_list[p].begin(), adj_list[p].end());
    }
    std::sort(adj_list[c].begin(), adj_list[c].end());
    
//    std::cout << "here5" << std::endl;
    for (unsigned int i = 0; i < adj_list[c].size(); i++){
        NodeID z  = adj_list[c][i];
        scratch[z] = false;
    }
//    std::cout << "here6" << std::endl;
    removeVertex(G, adj_list, a);
    removeVertex(G, adj_list, b);
//    std::cout << "here7" << std::endl;
    twin_fold.push_back(c);
//    std::cout << "here8" << std::endl;
}


void Reduction::removeTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c){
    
    std::vector<NodeID> new_clique1;
    new_clique1.push_back(v);
    new_clique1.push_back(a);
    new_clique1.push_back(b);
    addClique(G, adj_list, new_clique1);
    
    std::vector<NodeID> new_clique2;
    new_clique2.push_back(u);
    new_clique2.push_back(c);
    addClique(G, adj_list, new_clique2);
    
}


bool Reduction::removeTypeTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &a, NodeID &b, NodeID &c){
    
    for (unsigned int i = 0; i < adj_list[a].size(); i++){
        if (adj_list[a][i] == b){
            return true;
        }
        if (adj_list[a][i] == c){
            NodeID x = b;
            b = c;
            c = x;
            return true;
        }
    }
    for (unsigned int i = 0; i < adj_list[b].size(); i++){
        if (adj_list[b][i] == c){
            NodeID x = a;
            a = c;
            c = x;
            return true;
        }
    }
    return false;
}


bool Reduction::foundTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c){
    
//    std::cout << "vec: " <<  v << std::endl;
    for (unsigned int i = 0; i < adj_list[a].size(); i++){
        NodeID p = adj_list[a][i];
//        std::cout << "p: " <<  p << std::endl;
        if (p == v){
            continue;
        }
        if (adj_list[p].size() != 3){
            continue;
        }
        for (unsigned int j = 0; j < adj_list[b].size(); j++){
            NodeID q = adj_list[b][j];
//            std::cout << "q: " <<  q << std::endl;
            if (q > p) {
                break;
            }
            if (q == p){
                for (unsigned int k = 0; k < adj_list[c].size(); k++){
                    NodeID r = adj_list[c][k];
//                    std::cout << "r: " <<  r << std::endl;
                    if (r > p) {
                        break;
                    }
                    if (r == p){
                        u = p;
                        return true;
                    }
                }
                break;
            }
        }
    }
    return false;
}


void Reduction::organizeTwinElements(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &w, NodeID &y, NodeID &z, NodeID &a, NodeID &b, NodeID &c) {
    
    if (adj_list[w].size() < adj_list[y].size() && adj_list[w].size() < adj_list[z].size()){
        a = w;
        if (adj_list[y].size() < adj_list[z].size()) {
            b = y;
            c = z;
        }
        else {
            b = z;
            c = y;
        }
    }
    else if (adj_list[y].size() < adj_list[w].size() && adj_list[y].size() < adj_list[z].size()){
        a = y;
        if (adj_list[w].size() < adj_list[z].size()) {
            b = w;
            c = z;
        }
        else {
            b = z;
            c = w;
        }
    }
    else {
        a = z;
        if (adj_list[w].size() < adj_list[y].size()) {
            b = w;
            c = y;
        }
        else {
            b = y;
            c = w;
        }
    }
}


bool Reduction::checkValidTwinFold(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &a, NodeID &b, NodeID &c){
    
    NodeID w = adj_list[v][0];
    NodeID y = adj_list[v][1];
    NodeID z = adj_list[v][2];
//    std::cout << "heree" << std::endl;
    organizeTwinElements(G, adj_list, w, y, z, a, b, c);
    
    //        std::cout << w << ", " << y << ", " << z << std::endl;
    //        std::cout << a << ", " << b << ", " << c << std::endl;
//    std::cout << "hereee" << std::endl;
    if (foundTwin(G, adj_list, v, u, a, b, c)){
//        std::cout << "found valid" << std::endl;
//        std::cout << v << std::endl;
//        std::cout << w << ", " << y << ", " << z << std::endl;
//        std::cout << u << std::endl;
//        std::cout << adj_list[u][0] << ", " << adj_list[u][1] << ", " << adj_list[u][2] << std::endl;
        return true;
    }
    
    return false;
}


int Reduction::getTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    int num_folds = 0;
    forall_nodes(G, v) {
//        std::cout << v << std::endl;
        if (!node_status[v]){
//            std::cout << "seg" << std::endl;
            continue;
        }
        if (adj_list[v].size() == 3){
//            std::cout << "segg" << std::endl;
            NodeID u;
            NodeID a;
            NodeID b;
            NodeID c;
            //            std::cout << "found" << std::endl;
            if (checkValidTwinFold(G, adj_list, v, u, a, b, c)){
//                                std::cout << "valid twin" << std::endl;
//                std::cout << a << " " << b << " " << c;
                if (removeTypeTwin(G, adj_list, a, b, c)){
//                    std::cout << "remove" << std::endl;
//                    std::cout << a << " " << b << " " << c;
                    removeTwin(G, adj_list, v, u, a, b, c);
                }
                else {
                    foldTwinVertices(G, adj_list, v, u, a, b, c);
                    num_folds++;
//                std::cout << "seg?" << std::endl;
                    reductions.push_back(2);
                }
            }
        }
    } endfor
    return num_folds;
}


void Reduction::reduceTwin(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    bool folding = true;

    while (folding){
        int num_folds = getTwin(G, adj_list);
        if (num_folds == 0){
            folding = false;
        }
    }
//    getTwin(G, adj_list);
}


void Reduction::foldDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v){
    
    NodeID u = adj_list[v][0];
    NodeID w = adj_list[v][1];
    
    std::vector<NodeID> folded_nodes;
    folded_nodes.push_back(u);
    folded_nodes.push_back(v);
    folded_nodes.push_back(w);
    d2_folded_nodes.push_back(folded_nodes);
    
    d2_folded_neighborhood.push_back(adj_list[u]);
    
    removeVertex(G, adj_list, v);
    
    NodeID x = u;
    NodeID y = w;
    
    if (adj_list[u].size() < adj_list[w].size()){
        x = w;
        y = u;
    }
    
    for (unsigned int i = 0; i < adj_list[x].size(); i++){
        NodeID z  = adj_list[x][i];
        scratch[z] = true;
    }
    
    for (unsigned int i = 0; i < adj_list[y].size(); i++){
        NodeID z = adj_list[y][i];
        if (z == x){
            continue;
        }
        if (scratch[z]){
            continue;
        }
        adj_list[x].push_back(z);
        adj_list[z].push_back(x);
        std::sort(adj_list[z].begin(), adj_list[z].end());
    }
    std::sort(adj_list[x].begin(), adj_list[x].end());
    
    for (unsigned int i = 0; i < adj_list[x].size(); i++){
        NodeID z  = adj_list[x][i];
        scratch[z] = false;
    }
    
    removeVertex(G, adj_list, y);
    d2_fold.push_back(x);
}


bool Reduction::checkValidDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v){
    
    NodeID u = adj_list[v][0];
    NodeID w = adj_list[v][1];
    
    if (adj_list[u].size() == 2 && adj_list[w].size() == 2){
//        std::cout << adj_list[u].size() << adj_list[w].size() << std::endl;
        for (unsigned int i = 0; i < adj_list[u].size(); i++){
            if (adj_list[u][i] == w){
                return false;
            }
        }
    }
    
    for (unsigned int i = 0; i < adj_list[w].size(); i++){
        NodeID x = adj_list[w][i];
        scratch[x] = true;
    }
    for (unsigned int i = 0; i < adj_list[u].size(); i++){
        NodeID x = adj_list[u][i];
        scratch[x] = false;
    }
    
    for (unsigned int i = 0; i < adj_list[u].size(); i++){
        NodeID x = adj_list[u][i];
        
        if (!node_status[x]){
            continue;
        }
        if (x == v || x == w) {
            continue;
        }
        for (unsigned int j = 0; j < adj_list[x].size(); j++){
            NodeID y = adj_list[x][j];
            if (!node_status[y]){
                continue;
            }
            if (y == u || y == w){
                continue;
            }
            if (scratch[y]){
                for (unsigned int i = 0; i < adj_list[w].size(); i++){
                    NodeID x = adj_list[w][i];
                    scratch[x] = false;
                }
                return false;
            }
        }
    }
    for (unsigned int i = 0; i < adj_list[w].size(); i++){
        NodeID x = adj_list[w][i];
        scratch[x] = false;
    }
    
    return true;
}


int Reduction::getDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    int num_folds = 0;
    forall_nodes(G, v) {
        
        if (!node_status[v]){
            continue;
        }
        
        if (adj_list[v].size() == 2){
            // NodeID u = adj_list[v][0];
            // NodeID w = adj_list[v][1];
            if (checkValidDegreeTwo(G, adj_list, v)){
                foldDegreeTwo(G, adj_list, v);
                num_folds++;
                reductions.push_back(1);
//                std::cout << "valid: " << v << std::endl;
            }
        }
    } endfor
    return num_folds;
}


void Reduction::reduceDegreeTwo(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    bool folding = true;
    
    while (folding){
        int num_folds = getDegreeTwo(G, adj_list);
        if (num_folds == 0){
            folding = false;
        }
    }
}

//void followPath(graph_access &G, std::vector<bool> &node_status, std::vector<bool> &no\
//de_checked, NodeID &v){
//
//  node_checked[v] = true;
//
//  std::vector<NodeID> path_components;
//  path_components.push_back(v);
//
//  int path_size = 1;
//
//  EdgeID e = G.get_first_edge(v);
//
//  while (e != G.get_first_invalid_edge(v)) {
//  //forall_out_edges(G, e, v){
//    NodeID u = G.getEdgeTarget(e);
//
//    if (!node_status[u]) {
//    e++;
//    continue;
//    }
//
//    while (G.getNodeDegree(u) == 2 && !node_checked[u]){
//      node_checked[u] = true;
//      path_size++;
//
//      EdgeID e_u = G.get_first_edge(u);
//
//      if (!node_checked[G.getEdgeTarget(e_u)]){
//        u = G.getEdgeTarget(e_u);
//      }
//      else {
//        u = G.getEdgeTarget(e_u + 1);
//      }
//    }
//    e++;
//  }
//
//  std::cout << path_size << std::endl;
//}
//
//void getPath(graph_access &G, std::vector<bool> &node_status, std::vector<bool> &node_\
//checked, std::vector<NodeID> &path_vertices){
//
//  forall_nodes(G, v){
//
//    if (node_status[v] == false) {
//      continue;
//    }
//
//    if (G.getNodeDegree(v) == 2 && !node_checked[v]) {
//      std::cout << "Degree 2 vertex : " << v << std::endl;
//
//      followPath(G, node_status, node_checked, v);
//
//      path_vertices.push_back(v);
//    }
//  } endfor
//}


void Reduction::removeIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v){
    
    while (adj_list[v].size() > 0){
        NodeID u = adj_list[v][0];
        removeVertex(G, adj_list, u);
    }
    removeVertex(G, adj_list, v);
}


void Reduction::addIsolatedClique(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v){
    
    std::vector<NodeID> clique;
    clique.push_back(v);
    
    for (unsigned int i = 0; i < adj_list[v].size(); i++){
        NodeID u = adj_list[v][i];
        clique.push_back(u);
    }
    
    std::sort(clique.begin(), clique.end());
    clique_cover.push_back(clique);
    
    unsigned int cliqueID = clique_cover.size() - 1;
    for (unsigned int i = 0; i < clique.size(); i++){
        NodeID u = clique[i];
        node_clique[u] = cliqueID;
    }
}


bool Reduction::checkNeighbor(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u){
    
    EdgeID e_v = 0;
    EdgeID e_u = 0;
    
    bool not_found = false;
    
    while (e_v < adj_list[v].size()) {
        NodeID x = adj_list[v][e_v];
        NodeID w = adj_list[u][e_u];
        
        if (x == u) {
            e_v++;
            continue;
        }
        
        if (e_u == adj_list[u].size()){
            not_found = true;
            break;
        }
        
        if (x == w){
            e_v++;
            e_u++;
        }
        else if (x > w) {
            e_u++;
        }
        else {
            not_found = true;
            break;
        }
    }
    
    return not_found;
}


int Reduction::getIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    int num_isolated = 0;
    
    forall_nodes(G, v) {
        bool isolated = true;
        
        if (node_status[v] == false){
            continue;
        }
        
        for (unsigned int i = 0; i < adj_list[v].size(); i++){
            NodeID u = adj_list[v][i];
            
            bool not_found = checkNeighbor(G, adj_list, v, u);
            
            if (not_found){
                isolated = false;
                break;
            }
        }
        
        if (isolated) {
            num_isolated++;
            addIsolatedClique(G, adj_list, v);
            removeIsolated(G, adj_list, v);
            reductions.push_back(0);
        }
    } endfor
    
    return num_isolated;
}


void Reduction::reduceIsolated(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    int isolated_found = true;
    
    while (isolated_found){
        
        int num_isolated = getIsolated(G, adj_list);
        
        if (num_isolated == 0){
            isolated_found = false;
        }
    }
}


int Reduction::getGraphSize(graph_access &G){
    
    int in_nodes = 0;
    
    for (unsigned int i = 0; i < node_status.size(); i++){
        if (node_status[i]){
            in_nodes++;
        }
    }
    return in_nodes;
}


void Reduction::makeSubGraph(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, PartitionConfig &partition_config){

  makeNewAdj(G, adj_list);

  unsigned int n = 0;
  unsigned long e = 0;

  for (unsigned int i = 0; i < int_adj_list.size(); i++){
    for (unsigned int j = 0; j < int_adj_list[i].size(); j++){
      e++;
    }
    n++;
  }

  std::cout << n << " " << e / 2 << std::endl;
    
    cli *cli_instance;
    
    cli_instance = new cli(partition_config.seed);
    cli_instance->start_cli(int_adj_list, n, e);
    
    for (unsigned int i = 0; i < cli_instance->clique_cover.size(); i++){
        for (unsigned int j = 0; j < cli_instance->clique_cover[i].size(); j++){
            std::cout << cli_instance->clique_cover[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
    addCrownCliques(G, adj_list, cli_instance->clique_cover);
    
    delete(cli_instance);
    
    return;
  
  std::ofstream subGfile;
  subGfile.open("sub.graph");
  subGfile << n << " " << e / 2 << " 00\n";

  for (unsigned int i = 0; i < int_adj_list.size() - 1; i++){
    for (unsigned int j = 0; j < int_adj_list[i].size() - 1; j++){
      subGfile << int_adj_list[i][j] + 1 << " ";
    }
    subGfile << int_adj_list[i].back() + 1 << "\n";
  }
  for (unsigned int i = 0; i < int_adj_list.back().size() - 1; i++){
    subGfile << int_adj_list.back()[i] << " ";
  }
  subGfile << int_adj_list.back().back();  
  subGfile.close();
}


void Reduction::reduceGraph(graph_access &G, std::vector<std::vector<NodeID>> &adj_list, PartitionConfig &partition_config){
    
   unsigned int num_nodes = G.number_of_nodes();

   unsigned int old_G = 0;
   unsigned int new_G = num_nodes;

   while (new_G != old_G){
       old_G = new_G;
       // std::cout << "cycle" << std::endl;
       reduceIsolated(G, adj_list);
       reduceDegreeTwo(G, adj_list);
       reduceCrown(G, adj_list);
       reduceTwin(G, adj_list);
       reduceDominance(G, adj_list);
       //reduceCrown(G, adj_list);
       new_G = getGraphSize(G);
   }

     // reduceCrown(G, adj_list);
     // reduceIsolated(G, adj_list);
     
     //   reduceCrown(G, adj_list);
   makeNewAdj(G, adj_list);
   if (int_adj_list.size() > 0){
     std::cout << "remaining v" << std::endl;
     makeSubGraph(G, adj_list, partition_config);
   }
}

void makeAdjList(graph_access &G, std::vector<std::vector<NodeID>> &adj_list){
    
    forall_nodes(G, v) {
        std::vector<NodeID> neighbors;
        forall_out_edges(G, e, v){
            NodeID u = G.getEdgeTarget(e);
            neighbors.push_back(u);
        } endfor
        adj_list.push_back(neighbors);
    } endfor
}


void analyzeGraph(std::string &filename, graph_access &G, std::vector<std::vector<NodeID>> &adj_list, Reduction &R, timer &t){
    
    std::cout << filename << ", ";
    
    int num_nodes = G.number_of_nodes();
    std::cout << num_nodes << ", ";
    
    int num_edges = G.number_of_edges();
    std::cout << num_edges << ", ";
    
    int remaining_vertices = 0;
    for (unsigned int i = 0; i < R.node_status.size(); i++){
        if (R.node_status[i]){
            remaining_vertices++;
        }
    }
    
    std::cout << remaining_vertices << ", ";
    
    std::cout << t.elapsed() << ", ";
    
    std::cout << R.clique_cover.size() << ", ";
    
    std::vector<bool> in_clique(G.number_of_nodes(), false);
    for (unsigned int i = 0; i < R.clique_cover.size(); i++){
        std::vector<NodeID> clique = R.clique_cover[i];
        for (unsigned int j = 0; j < clique.size(); j++){
            NodeID x = clique[j];
            if (in_clique[x]){
                std::cout << "Duplicate: " << x << ", ";
            }
            else {
                in_clique[x] = true;
            }
        }
    }
    for (unsigned int i = 0; i < in_clique.size(); i++){
        if (in_clique[i] == false){
            std::cout << "Error" << std::endl;
            return;
        }
    }
    std::cout << "Complete" << std::endl;
}

int main(int argn, char **argv) {
    
    PartitionConfig partition_config;
    std::string graph_filename;

    bool is_graph_weighted = false;
    bool suppress_output   = false;
    bool recursive         = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if(ret_code) {
        return 0;
    }

    std::streambuf* backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if(suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }

    partition_config.LogDump(stdout);
    graph_access G;

    timer t;
    graph_io::readGraphWeighted(G, graph_filename);
    //  std::cout << "io time: " << t.elapsed()  << std::endl;

    //  unsigned int num_nodes = G.number_of_nodes();

    timer s;
    std::vector<std::vector<NodeID>> adj_list;
    makeAdjList(G, adj_list);
    //  std::vector<bool> node_status(num_nodes, true);

    Reduction R(G);
//      reduceGraph(G, adj_list, node_status);
    R.reduceGraph(G, adj_list, partition_config);
      R.unreduceGraph(G, adj_list);
    analyzeGraph(graph_filename, G, adj_list, R, s);
}


