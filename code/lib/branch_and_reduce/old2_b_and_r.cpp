#include <iostream>
#include <fstream>

#include "b_and_r.h"

branch_and_reduce::branch_and_reduce(graph_access &G, PartitionConfig &partition_config) {

  reduVCC = redu_vcc(G, partition_config);
  visited_nodes.assign(G.number_of_nodes(), false);

  mis = 0;
  num_reductions = 0;
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


void branch_and_reduce::branch( graph_access &G, unsigned int num_folded_cliques) {
  /* Simple Branch and Reduce -- No Pruning / Ordering of Eumnerated Cliques*/

  // return if past time limit
  if (t.elapsed() >= partition_config.solver_time_limit) return;

  // run exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  // number of fold cliques so far
  num_folded_cliques += R.num_fold_cliques;

  // kernel size is 0
  if (reduVCC.remaining_nodes == 0) {

    // compute current cover size
    unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

    // check if current cover is a better cover
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;

      // build cover
      reduVCC.build_cover(G);

      // unwind reductions
      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
          reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

      // check for valid cover
      reduVCC.validateCover(G);

    }

    // undo last reduction
    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }

  // get next node in kernel
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) {
    next_node++;
  }


  // enumerate all cliques for next_node
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

void branch_and_reduce::integrated_mis_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis) {
  /* Estimates kernel cover size using an integrated MIS from KaMIS.
     As vertices are removed from G, the MIS is updated. */

   // return if time limit is past
   if (t.elapsed() >= partition_config.solver_time_limit) return;

  // apply exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);

  // adjust number of folded cliques and kernel mis
  num_folded_cliques += R.num_fold_cliques;
  curr_mis -= (R.num_cliques + R.num_fold_cliques);

  // current cover number
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

  // check to see if kernel fully reduced
  if (reduVCC.remaining_nodes == 0) {

    // checks if current cover is better
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      // build cover
      reduVCC.build_cover(G);

      // unwind reductions
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
    }

    // undo last reductions
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

  // get next node in G'
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) next_node++;

  // estimate cover by adding current cover size and kernel mis
  unsigned int estimated_cover_size = curr_cover_size + curr_mis;
  // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;

  // prune if greater than current result
  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    // std::cout << "prune" << std::endl;
    // undo last set of reductions
    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }

  //enumerate all cliques of next_node
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {

    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);

    // recurr
    prune_branch(G, num_folded_cliques, curr_mis);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);
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
  // curr_mis -= (R.num_cliques + R.num_fold_cliques);
  curr_mis -= (R.num_cliques + R.num_fold_cliques);
  num_reductions += R.num_reductions;

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

      // reduVCC.validateCover(G);

      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
        reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

      reduVCC.validateCover(G);

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


void branch_and_reduce::lower_bound_branch( graph_access &G, unsigned int num_folded_cliques) {


  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  num_folded_cliques += R.num_fold_cliques;
  // curr_mis -= (R.num_cliques + R.num_fold_cliques);
  num_reductions += R.num_reductions;

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

      // reduVCC.validateCover(G);

      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
        reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

      reduVCC.validateCover(G);

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

  unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
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
    // curr_mis -= overlap;
    // std::cout << "branch" << std::endl;

    lower_bound_branch(G, num_folded_cliques);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    // curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC);
  reducer_stack.pop_back();
}


std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x) {

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(x);
  // std::cout << "complete enumerate" << std::endl;

  // sort enumerated cliques
  std::vector<unsigned int> curr_cliques_indices;
  std::vector<bool> curr_clique_is;
  std::vector<unsigned int> curr_clique_sizes;

  for (int i =0; i < curr_cliques.size(); i++) {

    curr_cliques_indices.push_back(i);
    curr_clique_sizes.push_back(curr_cliques[i].size());

    curr_clique_is.push_back(0);
    for (NodeID a : curr_cliques[i]) {
      if (reduVCC.node_mis[a]) {
        curr_clique_is.pop_back();
        curr_clique_is.push_back(1);
        break;
      }
    }
  }

  //order cliques
  std::sort(curr_cliques_indices.begin(), curr_cliques_indices.end(),
    [curr_clique_is, curr_clique_sizes](unsigned int i, unsigned int j) {


        if (curr_clique_is[i] == curr_clique_is[j]) {
          return curr_clique_sizes[i] < curr_clique_sizes[j];
        }
        return curr_clique_is[i] > curr_clique_is[j];
    });


    std::vector<std::vector<NodeID>> sorted_cliques;

    for (unsigned int index : curr_cliques_indices) {
      sorted_cliques.push_back(curr_cliques[index]);
    }

    return sorted_cliques;
}

void branch_and_reduce::sort_enumerate_branch( graph_access &G, unsigned int num_folded_cliques, PartitionConfig &partition_config, timer &t) {

  if (t.elapsed() >= partition_config.solver_time_limit) { return; }

  reducer R(G);
  R.exhaustive_reductions(G, reduVCC);
  reducer_stack.push_back(R);

  num_folded_cliques += R.num_fold_cliques;
  // curr_mis -= (R.num_cliques + R.num_fold_cliques);
  num_reductions += R.num_reductions;

  // increment next_node to find next node in the graph
  // check to see if we have made it through all nodes
  // if so, check to see if the cover is smaller, replace to min_cover

  if (reduVCC.remaining_nodes == 0) {
    unsigned int curr_cover_size = reduVCC.next_cliqueID;
    curr_cover_size += num_folded_cliques;
    if (curr_cover_size < reduVCC.clique_cover.size() || reduVCC.clique_cover.size() == 0) {

      // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);
      // R.unwindReductions(G, reduVCC);

      // reduVCC.validateCover(G);

      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
        reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

      reduVCC.validateCover(G);

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

  // std::cout << reduVCC.curr_mis << std::endl;
  unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
  // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;

  if (estimated_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    // std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC);
    reducer_stack.pop_back();
    return;
  }
  //enumerate all cliques
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;

  // branch on each clique
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G

    unsigned int overlap = overlapMIS(clique);

    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // curr_mis -= overlap;
    // std::cout << "branch" << std::endl;

    sort_enumerate_branch(G, num_folded_cliques, partition_config, t);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

    // curr_mis += overlap;
  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC);
  reducer_stack.pop_back();
}
