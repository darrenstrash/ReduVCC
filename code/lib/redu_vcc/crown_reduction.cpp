
#include <algorithm>
#include <iostream>
#include <fstream>

#include "crown_reduction.h"


void crown_reduction::generateAdjList(redu_vcc &reduVCC) {

  int_to_node_map.resize(reduVCC.num_nodes);
  node_to_int_map.resize(reduVCC.num_nodes);

  int node = 0;
  node_count = 0;

  for (NodeID v = 0; v < reduVCC.num_nodes; v++) {
    if (!reduVCC.node_status[v]) continue;
    int_to_node_map[node] = v;
    node_to_int_map[v] = node;
    node++; node_count++;
  }

  for (NodeID v = 0; v < reduVCC.num_nodes; v++) {
    std::vector<int> adj;
    if (!reduVCC.node_status[v]) continue;
    for (NodeID &u : reduVCC.adj_list[v]) {
      if (!reduVCC.node_status[u]) continue;
      int new_u = node_to_int_map[u];
      adj.push_back(new_u);
    }
    std::sort(adj.begin(), adj.end());
    adj_list.push_back(adj);
  }
}

void crown_reduction::addCliques(redu_vcc &reduVCC, std::vector<std::vector<int>> &cliques) {

  for (std::vector<int> &clique : cliques) {
    std::vector<NodeID> new_clique;
    for (int old_v : clique) {
      NodeID v = int_to_node_map[old_v];
      new_clique.push_back(v);
    }
    std::sort(new_clique.begin(), new_clique.end());
    crown_cliques.push_back(new_clique);
    reduVCC.addClique(new_clique);
    reduVCC.removeVertexSet(new_clique);
  }

}

void crown_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  type = "crown";

  generateAdjList(reduVCC);

  branch_and_reduce_algorithm b_and_r(adj_list, node_count);
  // reduVCC.buildKernel(G);
  //
  // branch_and_reduce_algorithm b_and_r(reduVCC.kernel_adj_list, reduVCC.remaining_nodes);
  // std::cout << "searched for crown" << std::endl;
  if (b_and_r.lpCrownReduction()){
    // std::cout << "crown found" << std::endl;
    addCliques(reduVCC, b_and_r.crown_cliques);
    unsigned int num_crown = crown_cliques.size();
    num_cliques += num_crown;

    // unsigned int curr_cliqueID = reduVCC.next_cliqueID;
    // reduVCC.addCrownCliques( crown_cliques, b_and_r.crown_cliques);
    // unsigned int num_crown = reduVCC.next_cliqueID - curr_cliqueID;
    // num_cliques += num_crown;
  }

  // std::cout << "here" << std::endl;
}

void crown_reduction::reduce( graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                           NodeID &node_v, NodeID &node_u ){
    reduce(G, reduVCC, node_v, node_u);

    for (std::vector<NodeID> clique : crown_cliques) {
      for (NodeID a : clique) queue->adjust_queue(reduVCC, a);;
    }
}

void crown_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){

  for (unsigned int i = crown_cliques.size(); i > 0; i--) {
    std::vector<NodeID> &clique = crown_cliques[i-1];
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);
  }
}
