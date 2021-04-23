
#include <iostream>
#include <fstream>

#include "iso_reduction.h"

bool iso_reduction::validNeighbor(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if u shares an edge with all w \in N(v)

    std::vector<std::vector<NodeID>>& adj_list = reduVCC.adj_list;

    // mark all vertices in N(v)
    unsigned int match_count = 0;
    for (NodeID w : adj_list[v]) {
      if (!reduVCC.node_status[w]) break;
      reduVCC.scratch1[w] = true;
      match_count++;
    }

    match_count--; // remove u from match count

    // check N(u) for matches
    for (NodeID w : adj_list[u]) {
      if (!reduVCC.node_status[w] || match_count == 0) break;
      if (reduVCC.scratch1[w]) match_count--;
    }

    for (NodeID w : adj_list[v]) reduVCC.scratch1[w] = false;

    if (match_count == 0) return true;
    return false;


    // unsigned int i = 0;
    // unsigned int j = 0;
    //
    // while (i < adj_list[v].size()) {
    //     NodeID x = adj_list[v][i];
    //     NodeID y = adj_list[u][j];
    //
    //     if (j == adj_list[u].size()){return false;}
    //
    //     if (!reduVCC.node_status[x]) { i++; continue; }
    //     if (!reduVCC.node_status[y]) { j++; continue; }
    //
    //     if (x == u) {i++; continue;}
    //
    //     if (x == y){i++; j++;}
    //     else if (x > y) {j++;}
    //     else {return false;}
    // }
    //
    // return true;
}

bool iso_reduction::validISO(redu_vcc &reduVCC, unsigned int &deg_limit, NodeID &v){
    // checks if v is an isolated vertex

    if (deg_limit != 0 && reduVCC.adj_size(v) > deg_limit) return false;
    // std::cout << reduVCC.adj_size(v) << std::endl;

    for (NodeID u : reduVCC.adj_list[v]) {
      if (!reduVCC.node_status[u]) break;
        // if (!reduVCC.node_status[u]) { continue; }

        bool valid_neighbor = iso_reduction::validNeighbor(reduVCC, v, u);
        if (!valid_neighbor) { return false; }
    }
    return true;

}

void iso_reduction::reduce(graph_access &G, redu_vcc &reduVCC,
                           NodeID &node_v, NodeID &node_u ){

  type = "iso";

  v = node_v;
  num_cliques++;

  deg = reduVCC.adj_size(v);

  // reduVCC.printNeighborhood(v);

  clique.push_back(v);
  for (NodeID u : reduVCC.adj_list[v]) {
    // if (!reduVCC.node_status[u]) { continue; };
    if (!reduVCC.node_status[u]) break;
    clique.push_back(u);
  };

  reduVCC.addClique(clique);
  reduVCC.removeVertexSet(clique);

}

void iso_reduction::reduce(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                           NodeID &node_v, NodeID &node_u ){
    reduce(G, reduVCC, node_v, node_u);
    for (NodeID a : clique) {
      if (a == v) continue;
      queue->adjust_queue(reduVCC, a);
    }

}

void iso_reduction::unreduce(graph_access &G, redu_vcc &reduVCC) {

  reduVCC.pop_clique(clique);
  reduVCC.addVertexSet(clique);
}

void iso_reduction::unfold(graph_access &G, redu_vcc &reduVCC) {
    // iso cliques do not need to be unfolded, G' is not impacted
}
