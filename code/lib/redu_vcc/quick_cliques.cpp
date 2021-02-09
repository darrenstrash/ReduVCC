
#include <algorithm>
#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "quick_cliques.h"


void quick_cliques::assignMaps(graph_access &G, redu_vcc &reduVCC) {

  old_to_new_map.clear();
  old_to_new_map.resize(G.number_of_nodes());
  new_to_old_map.clear();
  new_to_old_map.resize(G.number_of_nodes());

  int j = 0;
  for (unsigned int i = 0; i < G.number_of_nodes(); i++){
      if (!reduVCC.node_status[i]){
          continue;
      }
      old_to_new_map[i] = j;
      new_to_old_map[j] = i;
      j++;
  }
}

void quick_cliques::buildIntAdjList(graph_access &G, redu_vcc &reduVCC) {

  assignMaps(G);

  int_adj_list.clear();
  int_adj_list.resize(reduVCC.remaining_nodes);
  edges = 0;

  for (unsigned int i = 0; i < G.number_of_nodes(); i++) {
    if (!reduVCC.node_status[i]) { continue; }

    int new_v = old_to_new_map[i];

    std::list<int> adj;
    for (unsigned int j = 0; j < adj_list[i].size(); j++){
        NodeID u  = adj_list[i][j];
        int new_u = old_to_new_map[u];
        adj.push_back(new_u);
        kernel_edges++;
    }
    std::sort(adj.begin(), adj.end());
    int_adj_list[new_v] = adj;
  }
}
