
#include <algorithm>
#include <iostream>
#include <fstream>

#include "crown_reduction.h"



void crown_reduction::reduce( graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ){

  reduVCC.buildKernel(G);

  branch_and_reduce_algorithm b_and_r(reduVCC.kernel_adj_list, reduVCC.remaining_nodes);
  std::cout << "searched for crown" << std::endl;
  if (b_and_r.lpCrownReduction()){
    std::cout << "crown found" << std::endl;
    reduVCC.addCrownCliques( crown_cliques, b_and_r.crown_cliques);
  }
}

void crown_reduction::unreduce(graph_access &G, redu_vcc &reduVCC){

  for (unsigned int i = crown_cliques.size(); i > 0; i--) {
    std::vector<NodeID> &clique = crown_cliques[i-1];
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);
  }
}
