#include <iostream>
#include <fstream>

#include "b_and_r.h"

branch_and_reduce::branch_and_reduce(graph_access &G) {

  R = reducer(G);
  visited_nodes.assign(G.number_of_nodes(), false);

}

void branch_and_reduce::solveKernel(graph_access &G, PartitionConfig partition_config, timer &t) {
  R.solveKernel(G, partition_config, t);
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

    redu_vcc &reduVCC = R.reduVCC;

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

    pivot_enumerator(reduVCC, minimal_cliques, consider_nodes, curr_clique, excluded_nodes);

    // for (std::vector<NodeID> clique : minimal_cliques) R.reduVCC.printVectorSet(clique);
    return minimal_cliques;
}

void branch_and_reduce::pivot_enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
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


    pivot_enumerator(reduVCC, minimal_cliques, new_consider, curr_clique, new_excluded);

    curr_clique.pop_back();
    pivot_consider.erase(pivot_consider.begin());
    excluded_nodes.push_back(x);
  }


}

void branch_and_reduce::enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
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


    enumerator(reduVCC, minimal_cliques, new_consider, curr_clique, new_excluded);

    curr_clique.pop_back();
    consider_nodes.erase(consider_nodes.begin());
    excluded_nodes.push_back(x);
  }


}

unsigned int branch_and_reduce::remainingMIS(redu_vcc &reduVCC) {

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

unsigned int branch_and_reduce::exhaustive_reductions(graph_access &G, unsigned int &num_folded_cliques, unsigned int &curr_mis) {


  unsigned int all_reduced = 0;
  bool new_reduced = true;

  while (new_reduced) {
    new_reduced = false;

    std::vector<unsigned int> iso_stats = R.bruteISO(G);
    std::vector<unsigned int> d2_stats = R.bruteD2(G);
    std::vector<unsigned int> twin_stats = R.bruteTWIN(G);
    std::vector<unsigned int> dom_stats = R.bruteDOM(G);
    std::vector<unsigned int> crown_stats = R.bruteCROWN(G);

    unsigned int reduced = iso_stats[0] + d2_stats[0] + twin_stats[0] + dom_stats[0] + crown_stats[0];
    if (reduced) new_reduced = true;

    all_reduced += reduced;

    num_folded_cliques += (d2_stats[2] + twin_stats[2]);
    curr_mis -= (iso_stats[1] + d2_stats[1] + twin_stats[1] + dom_stats[1] + crown_stats[1] + num_folded_cliques);
  }

  return all_reduced;
}

void branch_and_reduce::branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node
                    ) {

  redu_vcc &reduVCC = R.reduVCC;

  unsigned int reduced = exhaustive_reductions(G, num_folded_cliques, curr_mis);

  NodeID next_node = curr_node;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover
  while (!reduVCC.node_status[next_node]) {
    if (next_node >= reduVCC.node_status.size()) {

      unsigned int curr_cover_size = reduVCC.next_cliqueID;
      curr_cover_size += num_folded_cliques;
      if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

        // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
        R.buildCover(G);
        R.unwindReductions(G);
      }

      R.undoReductions(G, reduced);
      return;
    }
    next_node++;
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

    branch(G, num_folded_cliques, curr_mis, next_node);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduced);
}


void branch_and_reduce::prune_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node
                    ) {

  redu_vcc &reduVCC = R.reduVCC;

  unsigned int reduced = exhaustive_reductions(G, num_folded_cliques, curr_mis);

  NodeID next_node = curr_node;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover
  while (!reduVCC.node_status[next_node]) {
    if (next_node >= reduVCC.node_status.size()) {

      unsigned int curr_cover_size = reduVCC.next_cliqueID;
      curr_cover_size += num_folded_cliques;
      if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

        // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
        R.buildCover(G);
        R.unwindReductions(G);
      }

      R.undoReductions(G, reduced);
      return;
    }
    next_node++;
  }


  // prune
  unsigned int curr_cover_size = reduVCC.next_cliqueID;
  curr_cover_size += num_folded_cliques;


  // how to prune:
  // if in a reduction -- remove 1 mis for each clique added
  // for a branch clique -- remove # independed verticies within the clique

  unsigned int estimated_cover_size = curr_cover_size + curr_mis;
  // std::cout << "est cover: " << estimated_cover_size << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    // std::cout << "prune" << std::endl;
    R.undoReductions(G, reduced);
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

    branch(G, num_folded_cliques, curr_mis, next_node);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduced);
}

void branch_and_reduce::small_deg_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node
                    ) {

  redu_vcc &reduVCC = R.reduVCC;

  unsigned int reduced = exhaustive_reductions(G, num_folded_cliques, curr_mis);
  // int reduced = 0;

  NodeID next_node = curr_node;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover

  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      R.buildCover(G);
      R.unwindReductions(G);
    }

    R.undoReductions(G, reduced);
    return;
  }

  unsigned int min_degree = 2;

  // while (!reduVCC.node_status[next_node]) {
  //   next_node++;
  // }

  while (true) {
    if (next_node >= reduVCC.node_status.size()) {
      min_degree++;
      next_node = curr_node;
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
  // std::cout << "est cover: " << estimated_cover_size << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    // std::cout << "prune" << std::endl;
    R.undoReductions(G, reduced);
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

    small_deg_branch(G, num_folded_cliques, curr_mis, next_node);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduced);
}

void branch_and_reduce::bounding_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node) {

  redu_vcc &reduVCC = R.reduVCC;

  unsigned int reduced = exhaustive_reductions(G, num_folded_cliques, curr_mis);
  // int reduced = 0;

  NodeID next_node = curr_node;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover
  std::cout << reduVCC.remaining_nodes << std::endl;
  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      R.buildCover(G);
      R.unwindReductions(G);
    }

    R.undoReductions(G, reduced);
    return;
  }

  unsigned int min_degree = 2;

  // while (!reduVCC.node_status[next_node]) {
  //   next_node++;
  // }

  while (true) {
    if (next_node >= reduVCC.node_status.size()) {
      min_degree++;
      next_node = curr_node;
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
  // std::cout << "est cover: " << estimated_cover_size << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
  // std::cout << "prune" << std::endl;
  R.undoReductions(G, reduced);
  return;
  }
  //enumerate all cliques
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;
  std::vector<unsigned int> curr_cliques_indices;
  std::vector<bool> curr_clique_is;
  std::vector<unsigned int> curr_clique_sizes;
  for (int i =0; i < curr_cliques.size(); i++) {
    curr_cliques_indices.push_back(i);
    curr_clique_sizes.push_back(curr_cliques[i].size());
    if (overlapMIS(curr_cliques[i])) {
      curr_clique_is.push_back(1);
    }
    else { curr_clique_is.push_back(0); }
  }
  //order cliques

  std::sort(curr_cliques_indices.begin(), curr_cliques_indices.end(),
    [curr_clique_is, curr_clique_sizes](unsigned int i, unsigned int j) {

        unsigned int s1 = curr_clique_is[i] * 1000 + curr_clique_sizes[i];
        unsigned int s2 = curr_clique_is[j] * 1000 + curr_clique_sizes[i];
        return s1 < s2;
    });

  for (unsigned int index : curr_cliques_indices) {
    std::cout << "mis: " << curr_clique_is[index] << ", size: " << curr_clique_sizes[index] << std::endl;
  }

  // branch on each clique
  // for (std::vector<NodeID> &clique : curr_cliques) {
  for (unsigned int index : curr_cliques_indices) {
  // add new clique and remove from G

  std::vector<NodeID> &clique = curr_cliques[index];

  // unsigned int overlap = overlapMIS(clique);
  unsigned int overlap = curr_clique_is[index];

  reduVCC.addClique(clique);
  reduVCC.removeVertexSet(clique);
  curr_mis -= overlap;
  // std::cout << "branch" << std::endl;

  bounding_branch(G, num_folded_cliques, curr_mis, next_node);

  // pop branched on clique
  reduVCC.pop_clique(clique);
  reduVCC.addVertexSet(clique);

  curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduced);
}



// experiments to do:
// no prune
// prune
// min degree prune

void branch_and_reduce::brute(graph_access &G) {

  redu_vcc &reduVCC = R.reduVCC;

  branch(G, 0, mis, 0);
}

void branch_and_reduce::prune(graph_access &G) {

  redu_vcc &reduVCC = R.reduVCC;

  prune_branch(G, 0, mis, 0);
}

void branch_and_reduce::min_degree_prune(graph_access &G) {

  redu_vcc &reduVCC = R.reduVCC;

  small_deg_branch(G, 0, mis, 0);
}

void branch_and_reduce::upper_bound_prune(graph_access &G, PartitionConfig &partition_config, timer &s){

  redu_vcc &reduVCC = R.reduVCC;
  R.solveKernel(G, partition_config, s);
  bounding_branch(G, 0, mis, 0);
}


void branch_and_reduce::analyzeGraph(std::string &filename, graph_access &G, timer &t) {
  R.analyzeGraph(filename, G, t);
}
