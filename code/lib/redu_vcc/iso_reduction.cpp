
#include <iostream>
#include <fstream>

#include "iso_reduction.h"

bool iso_reduction::validNeighbor(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if u shares an edge with all w \in N(v)

    std::vector<std::vector<NodeID>>& adj_list = reduVCC.adj_list;

    unsigned int i = 0;
    unsigned int j = 0;

    while (i < adj_list[v].size()) {
        NodeID x = adj_list[v][i];
        NodeID y = adj_list[u][j];

        if (x == u) {i++; continue;}

        if (j == adj_list[u].size()){return false;}

        if (x == y){i++; j++;}
        else if (x > y) {j++;}
        else {return false;}
    }

    return true;
}

bool iso_reduction::validISO(redu_vcc &reduVCC, NodeID &v){
    // checks if v is an isolated vertex

    for (NodeID u : reduVCC.adj_list[v]) {
        bool valid_neighbor = iso_reduction::validNeighbor(reduVCC, v, u);
        if (!valid_neighbor) { return false; }
    }
    return true;

}

void iso_reduction::reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  v = node_v;

  // reduVCC.printNeighborhood(v);

  std::vector<NodeID> clique;
  clique.push_back(v);
  for (NodeID u : reduVCC.adj_list[v]) {clique.push_back(u);};

  reduVCC.addClique(clique);
  reduVCC.removeVertexSet(clique);

}
