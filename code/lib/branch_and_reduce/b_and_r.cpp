
#include "b_and_r.h"

branch_and_reduce::branch_and_reduce(graph_access &G) {

  R.init(G);
  visited_nodes.assign(G.number_of_nodes(), false);
}

std::vector<std::vector<NodeID>> branch_and_reduce::enumerate(NodeID v) {

    std::vector<std::vector<NodeID>> minimal_cliques;

    std::vector<NodeID> consider_nodes;
    for (NodeID x : R.reduVCC.adj_list[v]) {
      if (!R.reduVCC.node_status[x]) { continue; }
      consider_nodes.push_back(x);
    }
    std::vector<NodeID> excluded_nodes {};
    std::vector<NodeID> curr_clique {v};

    enumerator(R.reduVCC, minimal_cliques, consider_nodes, curr_clique, excluded_nodes);

    return minimal_cliques;
}

void branch_and_reduce::enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
                                   std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes) {



  if (consider_nodes.empty() && excluded_nodes.empty()) {

    minimal_cliques.push_back(curr_clique);
  }

  while (!consider_nodes.empty()) {
    NodeID x = consider_nodes[0];

    for (NodeID y : reduVCC.adj_list[x]) {
      if (!reduVCC.node_status[y]) { continue; }
      reduVCC.scratch1[y] = true;
    }

    std::vector<NodeID> new_consider;
    for (NodeID y : consider_nodes) {
      if (reduVCC.scratch1[y]) { new_consider.push_back(y); }
    }

    std::vector<NodeID> new_excluded;
    for (NodeID y : excluded_nodes) {
      if (reduVCC.scratch1[y]) { new_excluded.push_back(y); }
    }

    for (NodeID y : reduVCC.adj_list[x]) {
      if (!reduVCC.node_status[y]) { continue; }
      reduVCC.scratch1[y] = false;
    }

    curr_clique.push_back(x);


    enumerator(reduVCC, minimal_cliques, new_consider, curr_clique, new_excluded);

    curr_clique.pop_back();
    consider_nodes.erase(consider_nodes.begin());
  }


}

void branch_and_reduce::branch( graph_access &G, unsigned int num_folded_cliques, NodeID curr_node
                    ) {

  redu_vcc &reduVCC = R.reduVCC;

  unsigned int i = R.bruteISO(G);
  unsigned int j = R.bruteD2(G);

  unsigned int reduced = i + j;
  num_folded_cliques += j;

  NodeID next_node = curr_node;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover
  while (!reduVCC.node_status[next_node]) {
    if (next_node >= reduVCC.node_status.size()) {

      unsigned int curr_cover_size = reduVCC.next_cliqueID;
      curr_cover_size += num_folded_cliques;
      if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

        std::cout << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
        R.buildCover(G);
        R.unwindReductions(G);
      }

      return;
    }
    next_node++;
  }

  //enumerate all cliques
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);

    reduVCC.removeVertexSet(clique);

    //branch on next clique
    branch(G, num_folded_cliques, next_node);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduced);
}


void branch_and_reduce::brute(graph_access &G) {

  redu_vcc &reduVCC = R.reduVCC;

  std::vector<std::vector<NodeID>> minimal_cover;

  branch(G, 0, 0);
  // for (std::vector<NodeID> &clique : minimal_cover) { reduVCC.printVectorSet(clique); }

}

void branch_and_reduce::analyzeGraph(std::string &filename, graph_access &G, timer &t) {
  R.analyzeGraph(filename, G, t);
}
