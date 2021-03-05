#include <iostream>
#include <fstream>

#include "b_and_r.h"

branch_and_reduce::branch_and_reduce(graph_access &G) {

  reduVCC = redu_vcc(G);
  visited_nodes.assign(G.number_of_nodes(), false);

  mis = 0;
}

void branch_and_reduce::getMIS(std::string file) {
  /* Generates node_mis mapping of minimum independent set from file */

  std::string line;

  std::ifstream mis_file (file);
  if (mis_file.is_open()) {
    while ( getline (mis_file, line)) {
      node_mis.push_back((int)line[0] - 48);
    }
    mis_file.close();
  }

  for (bool n : node_mis) if (n) { mis++; };
}

std::vector<std::vector<NodeID>> branch_and_reduce::enumerate(NodeID v) {
    /* Enumerates set of minimal cliques */

    // initialize set of minimal cliques
    std::vector<std::vector<NodeID>> minimal_cliques;

    // initial set of nodes to consider
    std::vector<NodeID> consider_nodes;
    for (NodeID x : reduVCC.adj_list[v]) {
      if (!reduVCC.node_status[x]) { continue; }
      consider_nodes.push_back(x);
    }

    // initial excluded notes and clique
    std::vector<NodeID> excluded_nodes {};
    std::vector<NodeID> curr_clique {v};

    pivot_enumerator( minimal_cliques, consider_nodes, curr_clique, excluded_nodes);

    // for (std::vector<NodeID> clique : minimal_cliques) R.reduVCC.printVectorSet(clique);
    return minimal_cliques;
}

void branch_and_reduce::pivot_enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                                   std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes) {


  if (consider_nodes.empty() && excluded_nodes.empty()) {

    minimal_cliques.push_back(curr_clique);
  }

  std::vector<NodeID> pivot_nodes = consider_nodes;
  NodeID pivot;
  std::vector<NodeID> N_pivot {};

  for (NodeID y : consider_nodes) reduVCC.scratch1[y] = true;
  for (NodeID y : excluded_nodes) if (!reduVCC.scratch1[y]) pivot_nodes.push_back(y);

  for (NodeID u : pivot_nodes) {
    unsigned int max = 0;
    for (NodeID y : reduVCC.adj_list[u]) {
      if (reduVCC.scratch1[y]) max++;
    }
    if (max > N_pivot.size()) {
      pivot = u;
      N_pivot = reduVCC.curr_adj_list(u);
      // std::cout << "pivot: " << u << " size: " << max << std::endl;
    }
  }

  for (NodeID y : N_pivot) reduVCC.scratch1[y] = false;


  // while (!consider_nodes.empty()) {
  std::vector<NodeID> pivot_consider {};
  for (NodeID x : consider_nodes) {
    if (reduVCC.scratch1[x]) pivot_consider.push_back(x);
  }

  for (NodeID y : consider_nodes) reduVCC.scratch1[y] = false;

  while (!pivot_consider.empty()) {


    // std::cout << consider_nodes.size() << std::endl;
    NodeID x = pivot_consider[0];

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


    pivot_enumerator(minimal_cliques, new_consider, curr_clique, new_excluded);

    curr_clique.pop_back();
    pivot_consider.erase(pivot_consider.begin());
    excluded_nodes.push_back(x);
  }


}

void branch_and_reduce::enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                                   std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes) {


  if (consider_nodes.empty() && excluded_nodes.empty()) {

    minimal_cliques.push_back(curr_clique);
  }

  while (!consider_nodes.empty()) {
    // std::cout << consider_nodes.size() << std::endl;
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


    enumerator(minimal_cliques, new_consider, curr_clique, new_excluded);

    curr_clique.pop_back();
    consider_nodes.erase(consider_nodes.begin());
    excluded_nodes.push_back(x);
  }


}

unsigned int branch_and_reduce::remainingMIS() {

  unsigned int kernelMIS = 0;

  for (unsigned int i = 0; i < node_mis.size(); i++) {
    if (!reduVCC.node_status[i]) continue;
    if (node_mis[i]) kernelMIS++;
  }

  return kernelMIS;
}

unsigned int branch_and_reduce::overlapMIS(std::vector<NodeID> &clique) {

  unsigned int overlap = 0;

  for (NodeID x : clique) if (node_mis[x]) overlap++;;

  return overlap;
}


// void branch_and_reduce::branch( graph_access &G, unsigned int num_fold_cliques, unsigned int curr_mis, NodeID curr_node
//                     ) {
//
//   reducer R(G);
//   R.exhaustive_reductions(G, reduVCC);
//
//   num_fold_cliques += R.num_fold_cliques;
//
//   NodeID next_node = curr_node;
//   // std::cout << next_node << std::endl;
//
//   // increment next_node to find next node in the graph
//   // check to see if we have made it through all nodes
//   // if so, check to see if the cover is smaller, replace to min_cover
//   while (!reduVCC.node_status[next_node]) {
//     if (next_node >= reduVCC.node_status.size()) {
//
//       unsigned int curr_cover_size = reduVCC.next_cliqueID;
//       curr_cover_size += num_fold_cliques;
//       if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {
//
//         std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//         reduVCC.build_cover(G);
//         R.unwindReductions(G, reduVCC);
//       }
//
//       R.undoReductions(G, reduVCC);
//       return;
//     }
//     next_node++;
//   }
//
//   //enumerate all cliques
//   // std::cout << "enumerate" << std::endl;
//   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
//   // std::cout << "complete enumerate" << std::endl;
//
//   // branch on each clique
//   for (std::vector<NodeID> &clique : curr_cliques) {
//     // add new clique and remove from G
//
//     unsigned int overlap = overlapMIS(clique);
//     // std::cout << "get overlap" << std::endl;
//
//     reduVCC.addClique(clique);
//     reduVCC.removeVertexSet(clique);
//     curr_mis -= overlap;
//     // std::cout << "branch" << std::endl;
//
//     branch(G, num_fold_cliques, curr_mis, next_node);
//
//     // pop branched on clique
//     reduVCC.pop_clique(clique);
//     reduVCC.addVertexSet(clique);
//
//     curr_mis += overlap;
//   }
//   // undo number of reductions from reduce
//   R.undoReductions(G, reduVCC);
// }


// void branch_and_reduce::prune_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node
//                     ) {
//
//   reducer R(G);
//   R.exhaustive_reductions(G, reduVCC);
//   reducer_stack.push_back(R);
//
//   num_folded_cliques += R.num_fold_cliques;
//   curr_mis -= R.num_cliques;
//
//   NodeID next_node = curr_node;
//
//   // increment next_node to find next node in the graph
//   // check to see if we have made it through all nodes
//   // if so, check to see if the cover is smaller, replace to min_cover
//
//
//   if (reduVCC.remaining_nodes == 0) {
//     unsigned int curr_cover_size = reduVCC.next_cliqueID;
//     curr_cover_size += num_folded_cliques;
//     if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {
//
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover(G);
//       // R.unwindReductions(G, reduVCC);
//
//
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) {
//           reducer_stack[i-1].unwindReductions(G, reduVCC);
//       }
//
//     }
//
//     R.undoReductions(G, reduVCC);
//     reducer_stack.pop_back();
//     return;
//   }
//
//   while (!reduVCC.node_status[next_node]) {
//     std::cout << next_node << std::endl;
//     next_node++;
//   }
//
//
//   // while (!reduVCC.node_status[next_node]) {
//   //   if (next_node >= reduVCC.node_status.size()) {
//   //
//   //     unsigned int curr_cover_size = reduVCC.next_cliqueID;
//   //     curr_cover_size += num_folded_cliques;
//   //     if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {
//   //
//   //       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   //       reduVCC.build_cover(G);
//   //       for (unsigned int i = reducer_stack.size(); i > 0; i--) {
//   //           reducer_stack[i-1].unwindReductions(G, reduVCC);
//   //       }
//   //     }
//   //
//   //     R.undoReductions(G, reduVCC);
//   //     reducer_stack.pop_back();
//   //     return;
//   //   }
//   //   next_node++;
//   // }
//
//
//   // prune
//   unsigned int curr_cover_size = reduVCC.next_cliqueID;
//   curr_cover_size += num_folded_cliques;
//
//
//   // how to prune:
//   // if in a reduction -- remove 1 mis for each clique added
//   // for a branch clique -- remove # independed verticies within the clique
//
//   unsigned int estimated_cover_size = curr_cover_size + curr_mis;
//   std::cout << "est cover: " << estimated_cover_size << ", "<<  reduVCC.clique_cover.size() << std::endl;
//
//   if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
//     // std::cout << "prune" << std::endl;
//     R.undoReductions(G, reduVCC);
//     reducer_stack.pop_back();
//     return;
//   }
//   //enumerate all cliques
//   // std::cout << "enumerate" << std::endl;
//   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
//   // std::cout << "complete enumerate" << std::endl;
//
//   // branch on each clique
//   for (std::vector<NodeID> &clique : curr_cliques) {
//     // add new clique and remove from G
//
//     unsigned int overlap = overlapMIS(clique);
//
//     reduVCC.addClique(clique);
//     reduVCC.removeVertexSet(clique);
//     curr_mis -= overlap;
//     // std::cout << "branch" << std::endl;
//
//     branch(G, num_folded_cliques, curr_mis, next_node);
//
//     // pop branched on clique
//     reduVCC.pop_clique(clique);
//     reduVCC.addVertexSet(clique);
//
//     curr_mis += overlap;
//   }
//   // undo number of reductions from reduce
//   R.undoReductions(G, reduVCC);
//   reducer_stack.pop_back();
// }

void branch_and_reduce::branch( graph_access &G, unsigned int num_folded_cliques) {


  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  num_folded_cliques += R.num_fold_cliques;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover

  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);
      // R.unwindReductions(G, reduVCC);


      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
          reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

    }

    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }

  NodeID next_node = 0;

  while (!reduVCC.node_status[next_node]) {
    next_node++;
  }


  //enumerate all cliques
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G

    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "branch" << std::endl;

    branch(G, num_folded_cliques);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC);
  reducer_stack.pop_back();
}

void branch_and_reduce::prune_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis
                    ) {


  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  num_folded_cliques += R.num_fold_cliques;
  curr_mis -= R.num_cliques;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover

  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);
      // R.unwindReductions(G, reduVCC);


      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
          reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

    }

    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }

  NodeID next_node = 0;

  while (!reduVCC.node_status[next_node]) {
    next_node++;
  }


  // prune
  unsigned int curr_cover_size = reduVCC.next_cliqueID;
  curr_cover_size += num_folded_cliques;


  // how to prune:
  // if in a reduction -- remove 1 mis for each clique added
  // for a branch clique -- remove # independed verticies within the clique

  unsigned int estimated_cover_size = curr_cover_size + curr_mis;
  std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }
  //enumerate all cliques
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G

    unsigned int overlap = overlapMIS(clique);

    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    curr_mis -= overlap;
    // std::cout << "branch" << std::endl;

    prune_branch(G, num_folded_cliques, curr_mis);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC);
  reducer_stack.pop_back();
}

void branch_and_reduce::small_deg_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis
                    ) {


  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  num_folded_cliques += R.num_fold_cliques;
  curr_mis -= R.num_cliques;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover

  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);
      // R.unwindReductions(G, reduVCC);


      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
          reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

    }

    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }

  unsigned int min_degree = 2;
  NodeID next_node = 0;

  // while (!reduVCC.node_status[next_node]) {
  //   next_node++;
  // }

  while (true) {
    if (next_node >= reduVCC.node_status.size()) {
      min_degree++;
      next_node = 0;
      continue;
    }
    if (!reduVCC.node_status[next_node]) {
      next_node++;
      continue;
    }

    if (reduVCC.adj_size(next_node) == min_degree) break;
    next_node++;
  }

  // prune
  unsigned int curr_cover_size = reduVCC.next_cliqueID;
  curr_cover_size += num_folded_cliques;


  // how to prune:
  // if in a reduction -- remove 1 mis for each clique added
  // for a branch clique -- remove # independed verticies within the clique

  unsigned int estimated_cover_size = curr_cover_size + curr_mis;
  std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }
  //enumerate all cliques
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G

    unsigned int overlap = overlapMIS(clique);

    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    curr_mis -= overlap;
    // std::cout << "branch" << std::endl;

    small_deg_branch(G, num_folded_cliques, curr_mis);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC);
  reducer_stack.pop_back();
}
//
//
// // experiments to do:
// // no prune
// // prune
// // min degree prune
//
// void branch_and_reduce::brute(graph_access &G) {
//
//   redu_vcc &reduVCC = R.reduVCC;
//
//   branch(G, 0, mis, 0);
// }
//
// void branch_and_reduce::prune(graph_access &G) {
//
//   redu_vcc &reduVCC = R.reduVCC;
//
//   prune_branch(G, 0, mis, 0);
// }
//
// void branch_and_reduce::min_degree_prune(graph_access &G) {
//
//   redu_vcc &reduVCC = R.reduVCC;
//
//   small_deg_branch(G, 0, mis, 0);
// }
//
// void branch_and_reduce::analyzeGraph(std::string &filename, graph_access &G, timer &t) {
//   R.analyzeGraph(filename, G, t);
// }
