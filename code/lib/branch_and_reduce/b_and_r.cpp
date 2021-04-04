#include <iostream>
#include <fstream>

#include "b_and_r.h"
#include "graph_io.h"
#include "mis/mis_config.h"
#include "mis/ils/ils.h"

branch_and_reduce::branch_and_reduce(graph_access &G) {

  reduVCC = redu_vcc(G);
  visited_nodes.assign(G.number_of_nodes(), false);

  num_reductions = 0;
  branch_count = 0;
}

branch_and_reduce::branch_and_reduce(graph_access &G, PartitionConfig &partition_config) {

  reduVCC = redu_vcc(G, partition_config);
  visited_nodes.assign(G.number_of_nodes(), false);

  num_reductions = 0;
  branch_count = 0;
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


void branch_and_reduce::brute_bandr( graph_access &G, unsigned int num_folded_cliques) {

  // perform exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
  // keep track of total folded cliques to determine current clique cover size
  num_folded_cliques += R.num_fold_cliques;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

  // get next node in kernel
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) next_node++;

  // enumerate all maximal cliques of next_node
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "branch" << std::endl;

    // branch
    branch_count++;
    brute_bandr(G, num_folded_cliques);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC); reducer_stack.pop_back();
}

void branch_and_reduce::reduMIS_bandr( graph_access &G, unsigned int num_folded_cliques) {

  // perform exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
  // keep track of total folded cliques to determine current clique cover size
  num_folded_cliques += R.num_fold_cliques;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

  // estimate cover size using partial cover size and MIS of kernel
  unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
  std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
  // prune branch if estimated cover is larger than current best
  if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
    std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }


  // get next node in kernel
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) next_node++;

  // enumerate all maximal cliques of next_node
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "branch" << std::endl;

    // branch
    branch_count++;
    reduMIS_bandr(G, num_folded_cliques);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC); reducer_stack.pop_back();
}

NodeID branch_and_reduce::min_deg_node(redu_vcc &reduVCC) {

  unsigned int min_degree = 0;
  NodeID next_node = 0;

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

    if (reduVCC.adj_size(next_node) == min_degree) return next_node;
    next_node++;
  }

}

void branch_and_reduce::small_degree_bandr( graph_access &G, unsigned int num_folded_cliques) {

  // perform exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
  // keep track of total folded cliques to determine current clique cover size
  num_folded_cliques += R.num_fold_cliques;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

  // estimate cover size using partial cover size and MIS of kernel
  unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
  std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
  // prune branch if estimated cover is larger than current best
  if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
    std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }


  // get next node in kernel with minimum degree
  NodeID next_node = min_deg_node(reduVCC);

  // enumerate all maximal cliques of next_node
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "branch" << std::endl;

    // branch
    branch_count++;
    small_degree_bandr(G, num_folded_cliques);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC); reducer_stack.pop_back();
}


std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x, std::vector<bool> &indset) {

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(x);
  //std::cout << "complete enumerate" << std::endl;

  // sort enumerated cliques
  std::vector<unsigned int> curr_cliques_indices;
  std::vector<bool> curr_clique_is;
  std::vector<unsigned int> curr_clique_sizes;

  for (int i =0; i < curr_cliques.size(); i++) {

    curr_cliques_indices.push_back(i);
    curr_clique_sizes.push_back(curr_cliques[i].size());

    curr_clique_is.push_back(0);
    for (NodeID a : curr_cliques[i]) {
      if (indset[a]) {
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
          return curr_clique_sizes[i] > curr_clique_sizes[j];
        }
        return curr_clique_is[i] < curr_clique_is[j];
    });


    std::vector<std::vector<NodeID>> sorted_cliques;

    for (unsigned int index : curr_cliques_indices) {
      sorted_cliques.push_back(curr_cliques[index]);
    }

    return sorted_cliques;
}

void branch_and_reduce::sort_enum_bandr( graph_access &G, unsigned int num_folded_cliques, PartitionConfig &partition_config, timer &t) {

  if (t.elapsed() > partition_config.solver_time_limit) return;

  // perform exhaustive reductions
  reducer R(G);
  R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
  // keep track of total folded cliques to determine current clique cover size
  num_folded_cliques += R.num_fold_cliques;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover(G);

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

  // estimate cover size using partial cover size and MIS of kernel
  unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
  // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
  // prune branch if estimated cover is larger than current best
  if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
    // std::cout << "prune" << std::endl;
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }


  // get next node in kernel with minimum degree
  NodeID next_node = min_deg_node(reduVCC);

  // enumerate all maximal cliques of next_node sorted by size and MIS
  // std::cout << "enumerate" << std::endl;
  std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "branch" << std::endl;

    // branch
    branch_count++;
    sort_enum_bandr(G, num_folded_cliques, partition_config, t);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC); reducer_stack.pop_back();
}

  void branch_and_reduce::chalupa_status_bandr( graph_access &G, unsigned int num_folded_cliques, PartitionConfig &partition_config, timer &t) {

    if (t.elapsed() > partition_config.solver_time_limit) return;

    // perform exhaustive reductions
    reducer R(G);
    R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
    // keep track of total folded cliques to determine current clique cover size
    num_folded_cliques += R.num_fold_cliques;

    // current size of parital clique cover
    unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;

    // check exit condition -- kernel is empty
    if (reduVCC.remaining_nodes == 0) {
      // check if we have a better solution
      if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
        // build current parital cover
        // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
        reduVCC.build_cover(G);

        // unwind reductions to get full cover
        for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
      }

      // undo branch's reductions and return
      R.undoReductions(G, reduVCC); reducer_stack.pop_back();
      return;
    }

    std::string graph_filename = "comparision";
    reduVCC.build_cover(G);
    reduVCC.solveKernel(G, partition_config, t);
    for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);

    // estimate cover size using partial cover size and MIS of kernel
    unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
    // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    // prune branch if estimated cover is larger than current best
    if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
      // std::cout << "prune" << std::endl;
      R.undoReductions(G, reduVCC); reducer_stack.pop_back();
      return;
    }


    // get next node in kernel with minimum degree
    NodeID next_node = min_deg_node(reduVCC);

    // enumerate all maximal cliques of next_node sorted by size and MIS
    // std::cout << "enumerate" << std::endl;
    std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
    // std::cout << "complete enumerate" << std::endl;
    // branch on each clique in enumerated set
    for (std::vector<NodeID> &clique : curr_cliques) {
      // add new clique and remove from G
      reduVCC.addClique(clique);
      reduVCC.removeVertexSet(clique);
      // std::cout << "branch" << std::endl;

      // branch
      branch_count++;
      chalupa_status_bandr(G, num_folded_cliques, partition_config, t);

      // pop branched on clique
      reduVCC.pop_clique(clique);
      reduVCC.addVertexSet(clique);

    }
    // undo number of reductions from reduce
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
  }

   void branch_and_reduce::generate_mis_bandr( graph_access &G, unsigned int num_folded_cliques, PartitionConfig &partition_config, timer &t) {
  
       if (t.elapsed() > partition_config.solver_time_limit) return;
  
       // perform exhaustive reductions
       reducer R(G);
       R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
       // keep track of total folded cliques to determine current clique cover size
       num_folded_cliques += R.num_fold_cliques;
  
       // current size of parital clique cover
       unsigned int curr_cover_size = reduVCC.next_cliqueID + num_folded_cliques;
  
       // check exit condition -- kernel is empty
       if (reduVCC.remaining_nodes == 0) {
         // check if we have a better solution
         if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
           // build current parital cover
           std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
           reduVCC.build_cover(G);
  
           // unwind reductions to get full cover
           for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
         }
  
         // undo branch's reductions and return
         R.undoReductions(G, reduVCC); reducer_stack.pop_back();
         return;
       }
 
       NodeID next_node;
       std::vector<std::vector<NodeID>> curr_cliques;

	 
       {
       // geneate MIS of kernel using ILS
       graph_access G_p;
       graph_io::readGraphKernel(G_p, reduVCC);
       MISConfig config;
       config.console_log = true;
       config.time_limit = 60;
       config.force_cand = 4;
       ils new_ils;
       new_ils.perform_ils(config, G_p, 1000);
  
  
       // estimate cover size using partial cover size and MIS of kernel
       unsigned int estimated_cover_size = curr_cover_size + new_ils.solution_size;
       // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
       // prune branch if estimated cover is larger than current best
       if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
         // std::cout << "prune" << std::endl;
         R.undoReductions(G, reduVCC); reducer_stack.pop_back();
         return;
       }
  
       // get next node in kernel with minimum degree
       next_node = min_deg_node(reduVCC);
  
       // enumerate all maximal cliques of next_node sorted by size and MIS
       // std::cout << "enumerate" << std::endl;
       curr_cliques = sorted_enumerate(next_node, new_ils.bool_solution);
       // std::cout << "complete enumerate" << std::endl;
       }

       // branch on each clique in enumerated set
       for (std::vector<NodeID> &clique : curr_cliques) {
         // add new clique and remove from G
         reduVCC.addClique(clique);
         reduVCC.removeVertexSet(clique);
         // std::cout << "branch" << std::endl;
  
         // branch
         branch_count++;
         generate_mis_bandr(G, num_folded_cliques, partition_config, t);
  
         // pop branched on clique
         reduVCC.pop_clique(clique);
         reduVCC.addVertexSet(clique);
  
       }
       // undo number of reductions from reduce
       R.undoReductions(G, reduVCC); reducer_stack.pop_back();
     }
