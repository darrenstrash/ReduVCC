#include <iostream>
#include <fstream>

#include "b_and_r.h"


void branch_and_reduce::construct_run(PartitionConfig &partition_config) {

  redu_type = "exhaustive";
  prune_type = "none";
  next_node_type = "none";
  enum_type = "none";

  if (partition_config.run_type == "brute") return;

  next_node_type = "small_deg";
  if (partition_config.run_type == "small_deg") return;
  enum_type = "sort_enum";
  if (partition_config.run_type == "sort_enum") return;
  prune_type = "ReduMIS";
  if (partition_config.run_type == "ReduMIS") return;
  redu_type = "cascading";
  if (partition_config.run_type == "cascading") return;
  prune_type = "KaMIS";
  if (partition_config.run_type == "KaMIS") return;
  prune_type = "SigMIS_nearlinear";
  if (partition_config.run_type == "SigMIS_nearlinear") return;
  prune_type = "SigMIS_linear";
  if (partition_config.run_type == "SigMIS_linear") return;
}

branch_and_reduce::branch_and_reduce(graph_access &G, redu_vcc &reduVCC, PartitionConfig &partition_config) {

  construct_run(partition_config);

  if (prune_type == "ReduMIS") reduVCC = redu_vcc(G, partition_config);
  else reduVCC = redu_vcc(G);

  branch_count = 0;
  prune_count = 0;
  decompose_count = 0;
  iso_degree.assign(reduVCC.num_nodes, 0);
  dom_degree.assign(reduVCC.num_nodes, 0);
  num_reductions = 0;
  num_attempts = 0;

  config.time_limit = 60;
  config.force_cand = 4;
}

branch_and_reduce::branch_and_reduce(redu_vcc &reduVCC, PartitionConfig &partition_config) {

  construct_run(partition_config);

  branch_count = 0;
  prune_count = 0;
  decompose_count = 0;
  iso_degree.assign(reduVCC.num_nodes, 0);
  dom_degree.assign(reduVCC.num_nodes, 0);
  num_reductions = 0;
  num_attempts = 0;

  config.time_limit = 60;
  config.force_cand = 4;
}

std::vector<std::vector<NodeID>> branch_and_reduce::enumerate(redu_vcc &reduVCC, NodeID v) {
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

// std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x, std::vector<bool> &indset) {
std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(redu_vcc &reduVCC, NodeID x) {

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(reduVCC, x);
  // std::cout << "complete enumerate " << std::endl;

  // sort enumerated cliques
  std::vector<unsigned int> curr_cliques_indices;
  // std::vector<bool> curr_clique_is;
  std::vector<unsigned int> curr_clique_sizes;

  for (int i =0; i < curr_cliques.size(); i++) {

    curr_cliques_indices.push_back(i);
    curr_clique_sizes.push_back(curr_cliques[i].size());

    // curr_clique_is.push_back(0);
    // for (NodeID a : curr_cliques[i]) {
    //   // std::cout << a << std::endl;
    //   if (indset[a]) {
    //     curr_clique_is.pop_back();
    //     curr_clique_is.push_back(1);
    //     break;
    //   }
    // }
  }

  // std::cout << "constructed" << std::endl;

  //order cliques
  std::sort(curr_cliques_indices.begin(), curr_cliques_indices.end(),
    // [curr_clique_is, curr_clique_sizes](unsigned int i, unsigned int j) {
    [curr_clique_sizes](unsigned int i, unsigned int j) {



        // if (curr_clique_is[i] == curr_clique_is[j]) {
          return curr_clique_sizes[i] > curr_clique_sizes[j];
        // }
        // return curr_clique_is[i] < curr_clique_is[j];
    });


    std::vector<std::vector<NodeID>> sorted_cliques;

    for (unsigned int index : curr_cliques_indices) {
      sorted_cliques.push_back(curr_cliques[index]);
    }

    return sorted_cliques;
}

void branch_and_reduce::reduce(redu_vcc &reduVCC, reducer &R, unsigned int &num_fold_cliques, vertex_queue *queue) {

    if (queue == nullptr || queue->empty()) {
      R.exhaustive_reductions(reduVCC, iso_degree, dom_degree);
    }
    else {
      R.cascading_reductions(reduVCC, queue, iso_degree, dom_degree);
    }
    reducer_stack.push_back(R);

    num_reductions += R.num_reductions;
    num_attempts += R.num_attempts;

    // keep track of total folded cliques to determine current clique cover size
    num_fold_cliques += R.num_fold_cliques;

}

bool branch_and_reduce::prune(redu_vcc &reduVCC, unsigned int &curr_cover_size) {

    unsigned int estimated_cover_size;

    if (prune_type == "none") {
      return false;
    }
    // else if (prune_type == "KaMIS") {
    //   // geneate MIS of kernel using ILS
    //   graph_access G_p;
    //   graph_io::readGraphKernel(G_p, reduVCC);
    //   // MISConfig config;
    //   // config.console_log = true;
    //   // config.time_limit = 60;
    //   // config.force_cand = 4;
    //   ils new_ils;
    //   new_ils.perform_ils(config, G_p, 1000);
    //
    //   estimated_cover_size = curr_cover_size + new_ils.solution_size;
    // }
    else if (prune_type == "SigMIS_linear") {
      // geneate MIS of kernel using Sigmod MIS
      Graph mis_G;
      mis_G.read_graph(reduVCC);
      unsigned int res_mis = mis_G.degree_two_kernal_and_remove_max_degree_without_contraction();

      estimated_cover_size = curr_cover_size + res_mis;
    }
    else if (prune_type == "SigMIS_nearlinear") {
      // geneate MIS of kernel using Sigmod MIS
      Graph mis_G;
      mis_G.read_graph(reduVCC);
      unsigned int res_mis = mis_G.degree_two_kernal_dominate_lp_and_remove_max_degree_without_contraction();

      estimated_cover_size = curr_cover_size + res_mis;
    }
    else { // prune_type == "ReduMIS"
      estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
    }

    // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    // prune branch if estimated cover is larger than current best
    if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
      // std::cout << "prune" << std::endl;
      prune_count++;
      return true;
    }

    return false;

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

std::vector<std::vector<NodeID>> branch_and_reduce::enum_vertex(redu_vcc &reduVCC, NodeID &v) {

  if (enum_type == "sort_enum") return sorted_enumerate(reduVCC, v);
  else return enumerate(reduVCC, v);
}

vertex_queue* branch_and_reduce::construct_queue(redu_vcc &reduVCC, std::vector<NodeID> &clique) {

  vertex_queue *new_queue = nullptr;
  if (redu_type == "exhaustive") return new_queue;

  new_queue = new vertex_queue(reduVCC);
  for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
  return new_queue;

}

NodeID branch_and_reduce::nextNode(redu_vcc &reduVCC){

  if (next_node_type == "small_deg") {
    return min_deg_node(reduVCC);
  }
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) next_node++;
  return next_node;

}

bool branch_and_reduce::decompose(redu_vcc &reduVCC, PartitionConfig &partition_config, timer &t,
                                  unsigned int curr_cover_size) {

  // std::cout << "remaining nodes: " << reduVCC.remaining_nodes << std::endl;
  std::vector<redu_vcc> children = reduVCC.decompose();

  if (children.empty()) return false;
  // std::cout << children.size() << std::endl;

  unsigned int cover_size = curr_cover_size;
  // std::cout << "curr cover: " << cover_size << std::endl;

  decompose_count += children.size();

  for (redu_vcc &child : children) {
    if (t.elapsed() > partition_config.solver_time_limit) return;

    // solve child
    // std::cout << "next child" << std::endl;
    branch_and_reduce B_child(child, partition_config);
    vertex_queue *queue = nullptr;
    if (partition_config.run_type == "cascading") queue = new vertex_queue(child);
    // child.printAdjList();
    B_child.bandr(child, 0, queue, partition_config, t);
    // std::cout << child.clique_cover.size() << std::endl;
    cover_size += child.clique_cover.size();

    // update branching counts
    branch_count += B_child.branch_count;
    prune_count += B_child.prune_count;
    // iso degree
    // dom_degree
    num_reductions += B_child.num_reductions;
    num_attempts += B_child.num_attempts;

  }

  if (t.elapsed() > partition_config.solver_time_limit) return;
  // std::cout << "cover size: " << cover_size << std::endl;

  // check clique cover size
  if (reduVCC.clique_cover.size() == 0 || cover_size < reduVCC.clique_cover.size()) {
    // for each child, add cliques to cover


    std::cout << "smaller cover: " << cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    reduVCC.build_cover();
    // std::cout << "builds?" << std::endl;
    for (redu_vcc &child : children) child.addCliquesToParent(reduVCC);

    // unwind reductions to get full cover
    for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
  }

  return true;
}

void branch_and_reduce::bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                               vertex_queue *queue,
                               PartitionConfig &partition_config, timer &t) {

  if (t.elapsed() > partition_config.solver_time_limit) return;

  reducer R(partition_config.iso_limit);
  reduce(reduVCC, R, num_fold_cliques, queue);
  delete queue;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
  std::cout << reduVCC.remaining_nodes << ", " << curr_cover_size << std::endl;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover();

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }


  if (prune(reduVCC, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }

  if (reduVCC.remaining_nodes > partition_config.decompose_limit && decompose(reduVCC, partition_config, t, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }

  // std::cout << "continues" << std::endl;
  // // estimate cover size using partial cover size and MIS of kernel
  // unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
  // // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
  // // prune branch if estimated cover is larger than current best
  // if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
  //   // std::cout << "prune" << std::endl;
  //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
  //   return;
  // }


  // get next node in kernel with minimum degree
  NodeID next_node = nextNode(reduVCC);

  // enumerate all maximal cliques of next_node sorted by size and MIS
  // std::cout << "enumerate" << std::endl;
  // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
  std::vector<std::vector<NodeID>> curr_cliques = enum_vertex(reduVCC, next_node);

  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "new queue" << std::endl;

    vertex_queue *new_queue = construct_queue(reduVCC, clique);
    // vertex_queue *new_queue = new vertex_queue(G);
    // for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);

    // std::cout << "branch" << std::endl;
    // branch
    branch_count++;
    bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(reduVCC); reducer_stack.pop_back();
}

// bool branch_and_reduce::edge_decompose(redu_vcc &reduVCC, PartitionConfig &partition_config, timer &t,
//                                   unsigned int curr_cover_size) {
//
//   // std::cout << "remaining nodes: " << reduVCC.remaining_nodes << std::endl;
//   std::vector<redu_vcc> children = reduVCC.decompose();
//
//   if (children.empty()) return false;
//   // std::cout << children.size() << std::endl;
//
//   unsigned int cover_size = curr_cover_size;
//   // std::cout << "curr cover: " << cover_size << std::endl;
//
//   decompose_count += children.size();
//
//   for (redu_vcc &child : children) {
//     if (t.elapsed() > partition_config.solver_time_limit) return;
//
//     // solve child
//     // std::cout << "next child" << std::endl;
//     branch_and_reduce B_child(child, partition_config);
//     vertex_queue *queue = nullptr;
//     if (partition_config.run_type == "cascading") queue = new vertex_queue(child);
//     // child.printAdjList();
//     B_child.bandr(child, 0, queue, partition_config, t);
//     // std::cout << child.clique_cover.size() << std::endl;
//     cover_size += child.clique_cover.size();
//
//     // update branching counts
//     branch_count += B_child.branch_count;
//     prune_count += B_child.prune_count;
//     // iso degree
//     // dom_degree
//     num_reductions += B_child.num_reductions;
//     num_attempts += B_child.num_attempts;
//
//   }
//
//   if (t.elapsed() > partition_config.solver_time_limit) return;
//   // std::cout << "cover size: " << cover_size << std::endl;
//
//   // check clique cover size
//   if (reduVCC.clique_cover.size() == 0 || cover_size < reduVCC.clique_cover.size()) {
//     // for each child, add cliques to cover
//
//
//     std::cout << "smaller cover: " << cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//     reduVCC.build_cover();
//     // std::cout << "builds?" << std::endl;
//     for (redu_vcc &child : children) child.addCliquesToParent(reduVCC);
//
//     // unwind reductions to get full cover
//     for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
//   }
//
//   return true;
// }

void branch_and_reduce::getEdge( redu_vcc &reduVCC, PartitionConfig &partition_config,
                                 NodeID &curr_node, NodeID &edge_node ) {

  while (true) {
   bool found_edge = false;

   if (!reduVCC.node_status[curr_node]) { curr_node++; continue; }
   for (NodeID u : reduVCC.adj_list[curr_node]) {
     if (reduVCC.node_status[u]) {
       edge_node = u;
       found_edge = true;
       break;
     }
   }
   if (found_edge) break;
  }

}

vertex_queue* branch_and_reduce::construct_queue_edge(redu_vcc &reduVCC, NodeID &v, NodeID &u) {

  vertex_queue *new_queue = nullptr;
  if (redu_type == "exhaustive") return new_queue;

  new_queue = new vertex_queue(reduVCC);
  new_queue->push(v);
  new_queue->adjust_queue(reduVCC, v);
  new_queue->adjust_queue(reduVCC, u);
  return new_queue;

}


void branch_and_reduce::branch_edge ( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                                    PartitionConfig &partition_config, timer &t,
                                    NodeID &curr_node, NodeID &edge_node ) {


  // save original curr_node neighborhood
  std::vector<NodeID> N_curr = reduVCC.adj_list[curr_node];

  // construct common neighborhood in curr_node, remove edge_node
  unsigned int i = 0;
  unsigned int j = 0;

  reduVCC.adj_list[curr_node].clear();
  while (i < N_curr.size() && j < reduVCC.adj_list[edge_node].size()) {
    NodeID v = N_curr[i];
    NodeID u = reduVCC.adj_list[edge_node][j];

    // if (!reduVCC.node_status[v]) { i++; continue; }
    // if (!reduVCC.node_status[u]) { j++; continue; }

    if (v < u) i++;
    else if (v > u) j++;
    else {
      reduVCC.adj_list[curr_node].push_back(v);
      j++; i++;
    }
  }

  reduVCC.removeVertex(edge_node);

  // mark nodes accordingly
  reduVCC.merge_node[curr_node] = true; // mark curr_node as a merge node
  reduVCC.merge_status[edge_node] = true;
  reduVCC.nodes_merged[curr_node].push_back(edge_node); // add edge_node to the set of merged nodes


  vertex_queue *new_queue = construct_queue_edge(reduVCC, curr_node, edge_node);

  edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);

  // undo branch
  reduVCC.nodes_merged[curr_node].pop_back(); // remove edge_node from set of merged nodes
  if (reduVCC.nodes_merged[curr_node].size() == 0) reduVCC.merge_node[curr_node] = false;
  reduVCC.merge_status[edge_node] = false;

  reduVCC.addVertex(edge_node); // add edge node back into graph
  reduVCC.adj_list[curr_node] = N_curr; // reset adj list
}

void branch_and_reduce::branch_nonedge ( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                                    PartitionConfig &partition_config, timer &t,
                                    NodeID &curr_node, NodeID &edge_node ) {

  // edge not in clique
  for (unsigned int i = 0; i < reduVCC.adj_list[curr_node].size(); i++) {
    if (reduVCC.adj_list[curr_node][i] == edge_node) {
      reduVCC.adj_list[curr_node].erase(reduVCC.adj_list[curr_node].begin() + i);
      break;
    }
  }
  for (unsigned int i = 0; i < reduVCC.adj_list[edge_node].size(); i++) {
    if (reduVCC.adj_list[edge_node][i] == curr_node) {
      reduVCC.adj_list[edge_node].erase(reduVCC.adj_list[edge_node].begin() + i);
      break;
    }
  }

  vertex_queue *new_queue = construct_queue_edge(reduVCC, curr_node, edge_node);
  // branch on no edge
  edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);

  reduVCC.adj_list[curr_node].push_back(edge_node);
  std::sort(reduVCC.adj_list[curr_node].begin(), reduVCC.adj_list[curr_node].end());
  reduVCC.adj_list[edge_node].push_back(curr_node);
  std::sort(reduVCC.adj_list[edge_node].begin(), reduVCC.adj_list[edge_node].end());



}


// void branch_and_reduce::edge_bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
//                                     vertex_queue *queue,
//                                     PartitionConfig &partition_config, timer &t,
//                                     NodeID curr_node
//                                   ) {
//
//   if (t.elapsed() > partition_config.solver_time_limit) return;
//
//   reducer R(partition_config.iso_limit);
//   reduce(reduVCC, R, num_fold_cliques, queue);
//   delete queue;
//
//   // current size of parital clique cover
//   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//   std::cout << reduVCC.remaining_nodes << ", " << curr_cover_size << std::endl;
//
//   // check exit condition -- kernel is empty
//   if (reduVCC.remaining_nodes == 0) {
//     // check if we have a better solution
//     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//       // build current parital cover
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover();
//
//       // unwind reductions to get full cover
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
//     }
//
//     // undo branch's reductions and return
//     R.undoReductions(reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//
//   if (prune(reduVCC, curr_cover_size)) {
//     R.undoReductions(reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//   // if (reduVCC.remaining_nodes > partition_config.decompose_limit && decompose(reduVCC, partition_config, t, curr_cover_size)) {
//   //   R.undoReductions(reduVCC); reducer_stack.pop_back();
//   //   return;
//   // }
//
//   // get next edge (adjusts curr_node if necessary)
//   NodeID edge_node; getEdge(reduVCC, partition_config, curr_node, edge_node);
//
//   // branch on edge
//   branch_edge(reduVCC, num_fold_cliques, partition_config, t, curr_node, edge_node);
//   // remove edge
//   branch_nonedge(reduVCC, num_fold_cliques, partition_config, t, curr_node, edge_node);
//
//
//   // undo number of reductions from reduce
//   R.undoReductions(reduVCC); reducer_stack.pop_back();
// }


void branch_and_reduce::edge_bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                               vertex_queue *queue,
                               PartitionConfig &partition_config, timer &t,
                               NodeID curr_node) {

  if (t.elapsed() > partition_config.solver_time_limit) return;


  reducer R(partition_config.iso_limit);
  reduce(reduVCC, R, num_fold_cliques, queue);
  delete queue;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
  // std::cout << reduVCC.remaining_nodes << ", " << curr_cover_size << std::endl;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    // std::cout << "cover size: " << curr_cover_size << std::endl;
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover();

      std::cout << "stack sizes: " << edge_stack.size() << ", " << reducer_stack.size() << std::endl;
      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) {
        reducer_stack[i-1].unwindReductions(reduVCC);
        if (i > 1) {
          std::vector<NodeID> &merged_edge = edge_stack[i-2];
          reduVCC.printVectorSet(merged_edge);
          unsigned int solveCliqueID = reduVCC.solve_node_clique[merged_edge[0]];
          reduVCC.solve_node_clique[merged_edge[1]] = solveCliqueID;
          reduVCC.clique_cover[solveCliqueID].push_back(merged_edge[1]);
          std::sort(reduVCC.clique_cover[solveCliqueID].begin(), reduVCC.clique_cover[solveCliqueID].end());
        }
      }

    }

    // undo branch's reductions and return
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }

  if (curr_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }


  if (prune(reduVCC, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return;
  }

  // if (reduVCC.remaining_nodes > partition_config.decompose_limit && decompose(reduVCC, partition_config, t, curr_cover_size)) {
  //   R.undoReductions(reduVCC); reducer_stack.pop_back();
  //   return;
  // }

  // if (decompose(reduVCC, partition_config, t, curr_cover_size)) {
  //   R.undoReductions(reduVCC); reducer_stack.pop_back();
  //   return;
  // }


  // std::cout << "reaches" << std::endl;
  // std::cout << std::endl << curr_cover_size << std::endl << std::endl;

  // checks if current node (node that is being merged is in graph)
  // NodeID edge_node;
  // while (true) {
  //   bool found_edge = false;
  //
  //   if (!reduVCC.node_status[curr_node]) { curr_node++; continue; }
  //   for (NodeID u : reduVCC.adj_list[curr_node]) {
  //     if (reduVCC.node_status[u]) {
  //       edge_node = u;
  //       found_edge = true;
  //       break;
  //     }
  //   }
  //   if (found_edge) break;
  // }

    // get next edge (adjusts curr_node if necessary)
    NodeID edge_node; getEdge(reduVCC, partition_config, curr_node, edge_node);

  // std::cout << "curr_node: " << curr_node << ", edge_node: " << edge_node << std::endl;


  // edge is in a clique

  // save original curr_node neighborhood
  std::vector<NodeID> N_curr = reduVCC.adj_list[curr_node];
  // construct common neighborhood in curr_node, remove edge_node
  unsigned int i = 0;
  unsigned int j = 0;

  reduVCC.adj_list[curr_node].clear();
  std::cout << reduVCC.adj_list[curr_node].size() << std::endl;
  while (i < N_curr.size() && j < reduVCC.adj_list[edge_node].size()) {
    NodeID v = N_curr[i];
    NodeID u = reduVCC.adj_list[edge_node][j];

  // if (!reduVCC.node_status[v]) { i++; continue; }
  // if (!reduVCC.node_status[u]) { j++; continue; }
  //
    if (v < u) {
      i++;
      // remove curr_node from neighborhood of v;
      for (unsigned int k = 0 ; k < reduVCC.adj_list[v].size(); k++) {
        if (reduVCC.adj_list[v][k] == curr_node) {
          reduVCC.adj_list[v].erase(reduVCC.adj_list[v].begin() + k);
          break;
        }
      }
    }

    else if (v > u) j++;
    else if (v == u){
      reduVCC.adj_list[curr_node].push_back(v);
      j++; i++;
    }
  }

  // std::cout << curr_node << ": ";
  // reduVCC.printVectorSet(N_curr);
  // std::cout << edge_node << ": ";
  // reduVCC.printVectorSet(reduVCC.adj_list[edge_node]);
  // reduVCC.printVectorSet(reduVCC.adj_list[curr_node]);
  // // std::cout << std::endl;
  //
  // for (NodeID a : reduVCC.adj_list[curr_node]) {
  //
  // }

  reduVCC.removeVertex(edge_node);

  // std::vector<NodeID> edge_clique {curr_node, edge_node};
  // add clique constructed by merge
  // if (reduVCC.merge_node[curr_node]) {  // if curr_node is already a merge
  //   unsigned int cliqueID = reduVCC.node_clique[curr_node]; // get ID
  //   reduVCC.node_clique[edge_node] = cliqueID;  // set edge_node ID
  // }
  // else reduVCC.addClique(edge_clique); // otherwise add new clique

  // reduVCC.addClique(edge_clique);

  reduVCC.merge_node[curr_node] = true; // mark curr_node as a merge node
  reduVCC.merge_status[edge_node] = true;
  reduVCC.nodes_merged[curr_node].push_back(edge_node); // add edge_node to the set of merged nodes


  std::vector<NodeID> edge {curr_node, edge_node};
  edge_stack.push_back(edge);
  // std::cout << "curr_node " << curr_node << ": ";
  // for (NodeID a : reduVCC.adj_list[curr_node]) {
  //   std::cout << a << ", ";
  // }
  // std::cout << std::endl;

  // branch on edge
  vertex_queue *new_queue = NULL;

  // std::cout << "edge: [" << curr_node << ", " << edge_node << "]" << std::endl;
  // std::cout << "edge add cover" << reduVCC.next_cliqueID +num_fold_cliques << std::endl;
  // std::cout << std::endl << std::endl;

  edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);

  // undo branch
  reduVCC.nodes_merged[curr_node].pop_back(); // remove edge_node from set of merged nodes
  // if (reduVCC.nodes_merged[curr_node].size() > 0) { reduVCC.node_clique[edge_node] = reduVCC.num_nodes; } // if merge has other nodes, remove edge_node from clique
  // else {  // only two in merged node
  //   reduVCC.pop_clique(edge_clique); // remove entire clique
  //   reduVCC.merge_node[curr_node] = false;  // curr_node is no longer merge node
  // }
  if (reduVCC.nodes_merged[curr_node].size() == 0) reduVCC.merge_node[curr_node] = false;
  reduVCC.merge_status[edge_node] = false;
  // reduVCC.pop_clique(edge_clique);

  // curr_node back to non common neighbors
  i = 0;
  j = 0;

  while (i < N_curr.size()) {
    NodeID v = N_curr[i];
    NodeID u = reduVCC.adj_list[curr_node][j];

    if (v == u) { // N[v] already contians curr_node
      i++; j++;
      continue;
    };
    // std::cout << "v: " << v << ", u:" << u << std::endl;
    reduVCC.insertVertex(v, curr_node);
    i++;
    // std::cout << i << std::endl;

  }


  // N_curr: 0 1 4 5 6 7 10
  // N[curr_node]: 1 5 7


  reduVCC.addVertex(edge_node); // add edge node back into graph
  reduVCC.adj_list[curr_node] = N_curr; // reset adj list

  edge_stack[-1].clear();


  // edge not in clique
  for (unsigned int i = 0; i < reduVCC.adj_list[curr_node].size(); i++) {
    if (reduVCC.adj_list[curr_node][i] == edge_node) {
      reduVCC.adj_list[curr_node].erase(reduVCC.adj_list[curr_node].begin() + i);
      break;
    }
  }
  for (unsigned int i = 0; i < reduVCC.adj_list[edge_node].size(); i++) {
    if (reduVCC.adj_list[edge_node][i] == curr_node) {
      reduVCC.adj_list[edge_node].erase(reduVCC.adj_list[edge_node].begin() + i);
      break;
    }
  }


  // branch on no edge
  edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);

  reduVCC.adj_list[curr_node].push_back(edge_node);
  std::sort(reduVCC.adj_list[curr_node].begin(), reduVCC.adj_list[curr_node].end());
  reduVCC.adj_list[edge_node].push_back(curr_node);
  std::sort(reduVCC.adj_list[edge_node].begin(), reduVCC.adj_list[edge_node].end());


  edge_stack.pop_back();




  // ---------------------- old branch vertex
  // // get next node in kernel with minimum degree
  // NodeID next_node = nextNode(reduVCC);
  //
  // // enumerate all maximal cliques of next_node sorted by size and MIS
  // // std::cout << "enumerate" << std::endl;
  // // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
  // std::vector<std::vector<NodeID>> curr_cliques = enum_vertex(reduVCC, next_node);
  //
  // // std::cout << "complete enumerate" << std::endl;
  // // branch on each clique in enumerated set
  // for (std::vector<NodeID> &clique : curr_cliques) {
  //   // add new clique and remove from G
  //   reduVCC.addClique(clique);
  //   reduVCC.removeVertexSet(clique);
  //   // std::cout << "new queue" << std::endl;
  //
  //   vertex_queue *new_queue = construct_queue(reduVCC, clique);
  //   // vertex_queue *new_queue = new vertex_queue(G);
  //   // for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
  //
  //   // std::cout << "branch" << std::endl;
  //   // branch
  //   branch_count++;
  //   bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t);
  //
  //   // pop branched on clique
  //   reduVCC.pop_clique(clique);
  //   reduVCC.addVertexSet(clique);
  //
  // }
    // ---------------------- old branch vertex



  // undo number of reductions from reduce
  R.undoReductions(reduVCC); reducer_stack.pop_back();
}



// void branch_and_reduce::edge_bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
//                                vertex_queue *queue,
//                                PartitionConfig &partition_config, timer &t,
//                                NodeID curr_node) {
//
//   // if (t.elapsed() > partition_config.solver_time_limit) return;
//   //
//
//   // // std::cout << "begin branch" << std::endl;
//   // reducer R(partition_config.iso_limit);
//   // reduce(reduVCC, R, num_fold_cliques, queue);
//   // // R.bruteISO(reduVCC, iso_degree);
//   // reducer_stack.push_back(R);
//   //
//   // num_reductions += R.num_reductions;
//   // num_attempts += R.num_attempts;
//   //
//   // // keep track of total folded cliques to determine current clique cover size
//   // num_fold_cliques += R.num_fold_cliques;
//   // delete queue;
//   //
//   // // current size of parital clique cover
//   // unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//   // std::cout << "nodes remaining: " << reduVCC.remaining_nodes << ", curr cover: " << curr_cover_size << std::endl;
//
//
//   reducer R(partition_config.iso_limit);
//   reduce(reduVCC, R, num_fold_cliques, queue);
//   delete queue;
//
//   // current size of parital clique cover
//   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//   std::cout << reduVCC.remaining_nodes << ", " << curr_cover_size << std::endl;
//
//   // check exit condition -- kernel is empty
//   if (reduVCC.remaining_nodes == 0) {
//     // check if we have a better solution
//     // std::cout << "cover size: " << curr_cover_size << std::endl;
//     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//       // build current parital cover
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover();
//
//       // unwind reductions to get full cover
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
//     }
//
//     // undo branch's reductions and return
//     R.undoReductions(reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//   if (curr_cover_size >= reduVCC.clique_cover.size() && reduVCC.clique_cover.size() != 0) {
//     R.undoReductions(reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//   //
//   // if (prune(reduVCC, curr_cover_size)) {
//   //   R.undoReductions(reduVCC); reducer_stack.pop_back();
//   //   return;
//   // }
//   //
//   // if (reduVCC.remaining_nodes > partition_config.decompose_limit && decompose(reduVCC, partition_config, t, curr_cover_size)) {
//   //   R.undoReductions(reduVCC); reducer_stack.pop_back();
//   //   return;
//   // }
//   //
//
//   // std::cout << "reaches" << std::endl;
//   // std::cout << std::endl << curr_cover_size << std::endl << std::endl;
//
//   // checks if current node (node that is being merged is in graph)
//   NodeID edge_node;
//   while (true) {
//     bool found_edge = false;
//
//     if (!reduVCC.node_status[curr_node]) { curr_node++; continue; }
//     for (NodeID u : reduVCC.adj_list[curr_node]) {
//       if (reduVCC.node_status[u]) {
//         edge_node = u;
//         found_edge = true;
//         break;
//       }
//     }
//     if (found_edge) break;
//   }
//
//   // std::cout << "curr_node: " << curr_node << ", edge_node: " << edge_node << std::endl;
//
//
//   // while (!reduVCC.node_status[curr_node]) curr_node++;
//   // NodeID edge_node;
//   // for (NodeID u : reduVCC.adj_list[curr_node]) {
//   //   if (reduVCC.node_status[u]) {
//   //     edge_node = u;
//   //     break;
//   //   }
//   // }
//
//   // std::cout << "curr_node " << curr_node << ": ";
//   // for (NodeID a : reduVCC.adj_list[curr_node]) {
//   //   if (reduVCC.node_status[a]) { std::cout << a << ", "; }
//   // }
//   // std::cout << std::endl;
//   // std::cout << "edge_node " << edge_node << ": ";
//   // for (NodeID a : reduVCC.adj_list[edge_node]) {
//   //   if (reduVCC.node_status[a]) { std::cout << a << ", "; }
//   // }
//   // std::cout << std::endl;
//
//   // edge is in a clique
//
//   // save original curr_node neighborhood
//   std::vector<NodeID> N_curr = reduVCC.adj_list[curr_node];
//   // construct common neighborhood in curr_node, remove edge_node
//   unsigned int i = 0;
//   unsigned int j = 0;
//
//   reduVCC.adj_list[curr_node].clear();
//   while (i < N_curr.size() && j < reduVCC.adj_list[edge_node].size()) {
//     NodeID v = N_curr[i];
//     NodeID u = reduVCC.adj_list[edge_node][j];
//
//   //   if (!reduVCC.node_status[u]) { j++; continue; }
//   //
//     if (v < u) i++;
//     else if (v > u) j++;
//     else {
//       reduVCC.adj_list[curr_node].push_back(v);
//       j++; i++;
//     }
//   }
//
//   reduVCC.removeVertex(edge_node);
//
//   // std::vector<NodeID> edge_clique {curr_node, edge_node};
//   // add clique constructed by merge
//   // if (reduVCC.merge_node[curr_node]) {  // if curr_node is already a merge
//   //   unsigned int cliqueID = reduVCC.node_clique[curr_node]; // get ID
//   //   reduVCC.node_clique[edge_node] = cliqueID;  // set edge_node ID
//   // }
//   // else reduVCC.addClique(edge_clique); // otherwise add new clique
//
//   // reduVCC.addClique(edge_clique);
//
//   reduVCC.merge_node[curr_node] = true; // mark curr_node as a merge node
//   reduVCC.merge_status[edge_node] = true;
//   reduVCC.nodes_merged[curr_node].push_back(edge_node); // add edge_node to the set of merged nodes
//
//   // std::cout << "curr_node " << curr_node << ": ";
//   // for (NodeID a : reduVCC.adj_list[curr_node]) {
//   //   std::cout << a << ", ";
//   // }
//   // std::cout << std::endl;
//
//   // branch on edge
//   vertex_queue *new_queue = NULL;
//
//   // std::cout << "edge: [" << curr_node << ", " << edge_node << "]" << std::endl;
//   // std::cout << "edge add cover" << reduVCC.next_cliqueID +num_fold_cliques << std::endl;
//   // std::cout << std::endl << std::endl;
//
//   edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);
//
//   // undo branch
//   reduVCC.nodes_merged[curr_node].pop_back(); // remove edge_node from set of merged nodes
//   // if (reduVCC.nodes_merged[curr_node].size() > 0) { reduVCC.node_clique[edge_node] = reduVCC.num_nodes; } // if merge has other nodes, remove edge_node from clique
//   // else {  // only two in merged node
//   //   reduVCC.pop_clique(edge_clique); // remove entire clique
//   //   reduVCC.merge_node[curr_node] = false;  // curr_node is no longer merge node
//   // }
//   if (reduVCC.nodes_merged[curr_node].size() == 0) reduVCC.merge_node[curr_node] = false;
//   reduVCC.merge_status[edge_node] = false;
//   // reduVCC.pop_clique(edge_clique);
//
//   reduVCC.addVertex(edge_node); // add edge node back into graph
//   reduVCC.adj_list[curr_node] = N_curr; // reset adj list
//
//
//   // edge not in clique
//   for (unsigned int i = 0; i < reduVCC.adj_list[curr_node].size(); i++) {
//     if (reduVCC.adj_list[curr_node][i] == edge_node) {
//       reduVCC.adj_list[curr_node].erase(reduVCC.adj_list[curr_node].begin() + i);
//       break;
//     }
//   }
//   for (unsigned int i = 0; i < reduVCC.adj_list[edge_node].size(); i++) {
//     if (reduVCC.adj_list[edge_node][i] == curr_node) {
//       reduVCC.adj_list[edge_node].erase(reduVCC.adj_list[edge_node].begin() + i);
//       break;
//     }
//   }
//
//
//   // branch on no edge
//   edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node);
//
//   reduVCC.adj_list[curr_node].push_back(edge_node);
//   std::sort(reduVCC.adj_list[curr_node].begin(), reduVCC.adj_list[curr_node].end());
//   reduVCC.adj_list[edge_node].push_back(curr_node);
//   std::sort(reduVCC.adj_list[edge_node].begin(), reduVCC.adj_list[edge_node].end());
//
//
//
//
//
//
//
//   // ---------------------- old branch vertex
//   // // get next node in kernel with minimum degree
//   // NodeID next_node = nextNode(reduVCC);
//   //
//   // // enumerate all maximal cliques of next_node sorted by size and MIS
//   // // std::cout << "enumerate" << std::endl;
//   // // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
//   // std::vector<std::vector<NodeID>> curr_cliques = enum_vertex(reduVCC, next_node);
//   //
//   // // std::cout << "complete enumerate" << std::endl;
//   // // branch on each clique in enumerated set
//   // for (std::vector<NodeID> &clique : curr_cliques) {
//   //   // add new clique and remove from G
//   //   reduVCC.addClique(clique);
//   //   reduVCC.removeVertexSet(clique);
//   //   // std::cout << "new queue" << std::endl;
//   //
//   //   vertex_queue *new_queue = construct_queue(reduVCC, clique);
//   //   // vertex_queue *new_queue = new vertex_queue(G);
//   //   // for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
//   //
//   //   // std::cout << "branch" << std::endl;
//   //   // branch
//   //   branch_count++;
//   //   bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t);
//   //
//   //   // pop branched on clique
//   //   reduVCC.pop_clique(clique);
//   //   reduVCC.addVertexSet(clique);
//   //
//   // }
//     // ---------------------- old branch vertex
//
//
//
//   // undo number of reductions from reduce
//   R.undoReductions(reduVCC); reducer_stack.pop_back();
// }




// bool branch_and_reduce::decompose() {
//
//   std::vector<redu_vcc> children = reduVCC
//
// }


// void branch_and_reduce::bandr( graph_access &G, unsigned int num_fold_cliques,
//                                vertex_queue *queue,
//                                PartitionConfig &partition_config, timer &t) {
//
//   if (t.elapsed() > partition_config.solver_time_limit) return;
//
//   reducer R(G);
//   if (queue->empty()) {
//     R.exhaustive_reductions(G, reduVCC, iso_degree, dom_degree);
//   }
//   else {
//     R.cascading_reductions(G, reduVCC, queue, iso_degree, dom_degree);
//   }
//   delete queue;
//   reducer_stack.push_back(R);
//
//   num_reductions += R.num_reductions;
//   num_attempts += R.num_attempts;
//
//   // keep track of total folded cliques to determine current clique cover size
//   num_fold_cliques += R.num_fold_cliques;
//
//   // current size of parital clique cover
//   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//
//   // check exit condition -- kernel is empty
//   if (reduVCC.remaining_nodes == 0) {
//     // check if we have a better solution
//     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//       // build current parital cover
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover(G);
//
//       // unwind reductions to get full cover
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
//     }
//
//     // undo branch's reductions and return
//     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//
//   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
//
// // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// // prune branch if estimated cover is larger than current best
// if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
//   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//   return;
// }
//
//   // // estimate cover size using partial cover size and MIS of kernel
//   // unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
//   // // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   // // prune branch if estimated cover is larger than current best
//   // if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
//   //   // std::cout << "prune" << std::endl;
//   //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//   //   return;
//   // }
//
//
//   // get next node in kernel with minimum degree
//   NodeID next_node = min_deg_node();
//
//   // enumerate all maximal cliques of next_node sorted by size and MIS
//   // std::cout << "enumerate" << std::endl;
//   // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
//   std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node);
//
//   // std::cout << "complete enumerate" << std::endl;
//   // branch on each clique in enumerated set
//   for (std::vector<NodeID> &clique : curr_cliques) {
//     // add new clique and remove from G
//     reduVCC.addClique(clique);
//     reduVCC.removeVertexSet(clique);
//     // std::cout << "new queue" << std::endl;
//
//     // vertex_queue *new_queue = construct_queue(G, clique);
//     vertex_queue *new_queue = new vertex_queue(G);
//     for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
//
//     // std::cout << "branch" << std::endl;
//     // branch
//     branch_count++;
//     bandr(G, num_fold_cliques, new_queue, partition_config, t);
//
//     // pop branched on clique
//     reduVCC.pop_clique(clique);
//     reduVCC.addVertexSet(clique);
//
//   }
//   // undo number of reductions from reduce
//   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// }


// // void branch_and_reduce::brute_bandr( graph_access &G, unsigned int num_fold_cliques) {
// //
// //   // perform exhaustive reductions
// //   reducer R(G);
// //   R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
// //   // keep track of total folded cliques to determine current clique cover size
// //   num_fold_cliques += R.num_fold_cliques;
// //
// //   // current size of parital clique cover
// //   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
// //
// //   // check exit condition -- kernel is empty
// //   if (reduVCC.remaining_nodes == 0) {
// //     // check if we have a better solution
// //     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
// //       // build current parital cover
// //       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //       reduVCC.build_cover(G);
// //
// //       // unwind reductions to get full cover
// //       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
// //     }
// //
// //     // undo branch's reductions and return
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //   // get next node in kernel
// //   NodeID next_node = 0;
// //   while (!reduVCC.node_status[next_node]) next_node++;
// //
// //   // enumerate all maximal cliques of next_node
// //   // std::cout << "enumerate" << std::endl;
// //   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
// //   // std::cout << "complete enumerate" << std::endl;
// //   // branch on each clique in enumerated set
// //   for (std::vector<NodeID> &clique : curr_cliques) {
// //     // add new clique and remove from G
// //     reduVCC.addClique(clique);
// //     reduVCC.removeVertexSet(clique);
// //     // std::cout << "branch" << std::endl;
// //
// //     // branch
// //     branch_count++;
// //     brute_bandr(G, num_fold_cliques);
// //
// //     // pop branched on clique
// //     reduVCC.pop_clique(clique);
// //     reduVCC.addVertexSet(clique);
// //
// //   }
// //   // undo number of reductions from reduce
// //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// // }
//
// // void branch_and_reduce::reduMIS_bandr( graph_access &G, unsigned int num_fold_cliques) {
// //
// //   // perform exhaustive reductions
// //   reducer R(G);
// //   R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
// //   // keep track of total folded cliques to determine current clique cover size
// //   num_fold_cliques += R.num_fold_cliques;
// //
// //   // current size of parital clique cover
// //   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
// //
// //   // check exit condition -- kernel is empty
// //   if (reduVCC.remaining_nodes == 0) {
// //     // check if we have a better solution
// //     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
// //       // build current parital cover
// //       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //       reduVCC.build_cover(G);
// //
// //       // unwind reductions to get full cover
// //       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
// //     }
// //
// //     // undo branch's reductions and return
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //   // estimate cover size using partial cover size and MIS of kernel
// //   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
// //   std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //   // prune branch if estimated cover is larger than current best
// //   if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
// //     std::cout << "prune" << std::endl;
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //
// //   // get next node in kernel
// //   NodeID next_node = 0;
// //   while (!reduVCC.node_status[next_node]) next_node++;
// //
// //   // enumerate all maximal cliques of next_node
// //   // std::cout << "enumerate" << std::endl;
// //   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
// //   // std::cout << "complete enumerate" << std::endl;
// //   // branch on each clique in enumerated set
// //   for (std::vector<NodeID> &clique : curr_cliques) {
// //     // add new clique and remove from G
// //     reduVCC.addClique(clique);
// //     reduVCC.removeVertexSet(clique);
// //     // std::cout << "branch" << std::endl;
// //
// //     // branch
// //     branch_count++;
// //     reduMIS_bandr(G, num_fold_cliques);
// //
// //     // pop branched on clique
// //     reduVCC.pop_clique(clique);
// //     reduVCC.addVertexSet(clique);
// //
// //   }
// //   // undo number of reductions from reduce
// //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// // }
//
// NodeID branch_and_reduce::min_deg_node(redu_vcc &reduVCC) {
//
//   unsigned int min_degree = 0;
//   NodeID next_node = 0;
//
//   while (true) {
//     if (next_node >= reduVCC.node_status.size()) {
//       min_degree++;
//       next_node = 0;
//       continue;
//     }
//     if (!reduVCC.node_status[next_node]) {
//       next_node++;
//       continue;
//     }
//
//     if (reduVCC.adj_size(next_node) == min_degree) return next_node;
//     next_node++;
//   }
//
// }
//
// void branch_and_reduce::small_degree_bandr( graph_access &G, unsigned int num_fold_cliques) {
//
//   // perform exhaustive reductions
//   reducer R(G);
//   R.exhaustive_reductions(G, reduVCC, iso_degree, dom_degree); reducer_stack.push_back(R);
//   red_perf += R.num_reductions;
//   red_tried += R.num_attempts;
//   // keep track of total folded cliques to determine current clique cover size
//   num_fold_cliques += R.num_fold_cliques;
//
//   // current size of parital clique cover
//   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//
//   // check exit condition -- kernel is empty
//   if (reduVCC.remaining_nodes == 0) {
//     // check if we have a better solution
//     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//       // build current parital cover
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover(G);
//
//       // unwind reductions to get full cover
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
//     }
//
//     // undo branch's reductions and return
//     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//   // estimate cover size using partial cover size and MIS of kernel
//   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
//   // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   // prune branch if estimated cover is larger than current best
//   if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
//     // std::cout << "prune" << std::endl;
//     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//
//   // get next node in kernel with minimum degree
//   NodeID next_node = min_deg_node(reduVCC);
//
//   // enumerate all maximal cliques of next_node
//   // std::cout << "enumerate" << std::endl;
//   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
//   // std::cout << "complete enumerate" << std::endl;
//   // branch on each clique in enumerated set
//   for (std::vector<NodeID> &clique : curr_cliques) {
//     // add new clique and remove from G
//     reduVCC.addClique(clique);
//     reduVCC.removeVertexSet(clique);
//     // std::cout << "branch" << std::endl;
//
//     // branch
//     branch_count++;
//     small_degree_bandr(G, num_fold_cliques);
//
//     // pop branched on clique
//     reduVCC.pop_clique(clique);
//     reduVCC.addVertexSet(clique);
//
//   }
//   // undo number of reductions from reduce
//   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// }
//
//
// std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x, std::vector<bool> &indset) {
//
//   std::vector<std::vector<NodeID>> curr_cliques = enumerate(x);
//   // std::cout << "complete enumerate " << std::endl;
//
//   // sort enumerated cliques
//   std::vector<unsigned int> curr_cliques_indices;
//   std::vector<bool> curr_clique_is;
//   std::vector<unsigned int> curr_clique_sizes;
//
//   for (int i =0; i < curr_cliques.size(); i++) {
//
//     curr_cliques_indices.push_back(i);
//     curr_clique_sizes.push_back(curr_cliques[i].size());
//
//     curr_clique_is.push_back(0);
//     for (NodeID a : curr_cliques[i]) {
//       // std::cout << a << std::endl;
//       if (indset[a]) {
//         curr_clique_is.pop_back();
//         curr_clique_is.push_back(1);
//         break;
//       }
//     }
//   }
//
//   // std::cout << "constructed" << std::endl;
//
//   //order cliques
//   std::sort(curr_cliques_indices.begin(), curr_cliques_indices.end(),
//     [curr_clique_is, curr_clique_sizes](unsigned int i, unsigned int j) {
//
//
//         if (curr_clique_is[i] == curr_clique_is[j]) {
//           return curr_clique_sizes[i] > curr_clique_sizes[j];
//         }
//         return curr_clique_is[i] < curr_clique_is[j];
//     });
//
//
//     std::vector<std::vector<NodeID>> sorted_cliques;
//
//     for (unsigned int index : curr_cliques_indices) {
//       sorted_cliques.push_back(curr_cliques[index]);
//     }
//
//     return sorted_cliques;
// }
//
// // void branch_and_reduce::sort_enum_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t) {
// //
// //   if (t.elapsed() > partition_config.solver_time_limit) return;
// //
// //   // perform exhaustive reductions
// //   reducer R(G);
// //   R.exhaustive_reductions(G, reduVCC, iso_degree, dom_degree); reducer_stack.push_back(R);
// //   // keep track of total folded cliques to determine current clique cover size
// //   num_fold_cliques += R.num_fold_cliques;
// //
// //   // current size of parital clique cover
// //   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
// //
// //   // check exit condition -- kernel is empty
// //   if (reduVCC.remaining_nodes == 0) {
// //     // check if we have a better solution
// //     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
// //       // build current parital cover
// //       // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //       reduVCC.build_cover(G);
// //
// //       // unwind reductions to get full cover
// //       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
// //     }
// //
// //     // undo branch's reductions and return
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //   // estimate cover size using partial cover size and MIS of kernel
// //   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
// //   // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //   // prune branch if estimated cover is larger than current best
// //   if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
// //     // std::cout << "prune" << std::endl;
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //
// //   // get next node in kernel with minimum degree
// //   NodeID next_node = min_deg_node(reduVCC);
// //
// //   // enumerate all maximal cliques of next_node sorted by size and MIS
// //   // std::cout << "enumerate" << std::endl;
// //   std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
// //   // std::cout << "complete enumerate" << std::endl;
// //   // branch on each clique in enumerated set
// //   for (std::vector<NodeID> &clique : curr_cliques) {
// //     // add new clique and remove from G
// //     reduVCC.addClique(clique);
// //     reduVCC.removeVertexSet(clique);
// //     // std::cout << "branch" << std::endl;
// //
// //     // branch
// //     branch_count++;
// //     sort_enum_bandr(G, num_fold_cliques, partition_config, t);
// //
// //     // pop branched on clique
// //     reduVCC.pop_clique(clique);
// //     reduVCC.addVertexSet(clique);
// //
// //   }
// //   // undo number of reductions from reduce
// //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// // }
//
// // void branch_and_reduce::chalupa_status_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t) {
// //
// //   if (t.elapsed() > partition_config.solver_time_limit) return;
// //
// //   // perform exhaustive reductions
// //   reducer R(G);
// //   R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
// //   // keep track of total folded cliques to determine current clique cover size
// //   num_fold_cliques += R.num_fold_cliques;
// //
// //   // current size of parital clique cover
// //   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
// //
// //   // check exit condition -- kernel is empty
// //   if (reduVCC.remaining_nodes == 0) {
// //     // check if we have a better solution
// //     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
// //       // build current parital cover
// //       // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //       reduVCC.build_cover(G);
// //
// //       // unwind reductions to get full cover
// //       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
// //     }
// //
// //     // undo branch's reductions and return
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //   std::string graph_filename = "comparision";
// //   reduVCC.build_cover(G);
// //   reduVCC.solveKernel(G, partition_config, t);
// //   for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
// //
// //   // estimate cover size using partial cover size and MIS of kernel
// //   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
// //   // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
// //   // prune branch if estimated cover is larger than current best
// //   if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
// //     // std::cout << "prune" << std::endl;
// //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// //     return;
// //   }
// //
// //
// //   // get next node in kernel with minimum degree
// //   NodeID next_node = min_deg_node(reduVCC);
// //
// //   // enumerate all maximal cliques of next_node sorted by size and MIS
// //   // std::cout << "enumerate" << std::endl;
// //   std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
// //   // std::cout << "complete enumerate" << std::endl;
// //   // branch on each clique in enumerated set
// //   for (std::vector<NodeID> &clique : curr_cliques) {
// //     // add new clique and remove from G
// //     reduVCC.addClique(clique);
// //     reduVCC.removeVertexSet(clique);
// //     // std::cout << "branch" << std::endl;
// //
// //     // branch
// //     branch_count++;
// //     chalupa_status_bandr(G, num_fold_cliques, partition_config, t);
// //
// //     // pop branched on clique
// //     reduVCC.pop_clique(clique);
// //     reduVCC.addVertexSet(clique);
// //
// //   }
// //   // undo number of reductions from reduce
// //   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// // }
//
// void branch_and_reduce::cascading_red_bandr( graph_access &G, unsigned int num_fold_cliques,
//                                              vertex_queue *queue,
//                                              PartitionConfig &partition_config, timer &t) {
//
//   if (t.elapsed() > partition_config.solver_time_limit) return;
//
//   // perform exhaustive reductions
//   // std::cout << "pre reductions reduMIS " << reduVCC.curr_mis << std::endl;
//   reducer R(G);
//   if (queue->empty()) {
//     R.exhaustive_reductions(G, reduVCC, iso_degree, dom_degree);
//   } else {
//     R.cascading_reductions(G, reduVCC, queue, iso_degree, dom_degree);
//   }
//   reducer_stack.push_back(R);
//   red_perf += R.num_reductions;
//   red_tried += R.num_attempts;
//   // R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
//   // keep track of total folded cliques to determine current clique cover size
//   num_fold_cliques += R.num_fold_cliques;
//   delete queue;
//
//   // current size of parital clique cover
//   unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//   // std::cout << "curr cover size: " << curr_cover_size << std::endl;
//   // std::cout << "post reductions reduMIS: " << reduVCC.curr_mis << std::endl;
//   // std::cout << "remaining nodes: " << reduVCC.remaining_nodes << std::endl;
//
//   // check exit condition -- kernel is empty
//   if (reduVCC.remaining_nodes == 0) {
//     // check if we have a better solution
//     if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//       // build current parital cover
//       std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//       reduVCC.build_cover(G);
//
//       // unwind reductions to get full cover
//       for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
//     }
//
//     // undo branch's reductions and return
//     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//   // estimate cover size using partial cover size and MIS of kernel
//   unsigned int estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
//   // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   // prune branch if estimated cover is larger than current best
//   if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
//     // std::cout << "prune" << std::endl;
//     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//     return;
//   }
//
//
//   // get next node in kernel with minimum degree
//   NodeID next_node = min_deg_node(reduVCC);
//
//   // enumerate all maximal cliques of next_node sorted by size and MIS
//   // std::cout << "enumerate" << std::endl;
//   // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
//   std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);
//
//   // std::cout << "complete enumerate" << std::endl;
//   // branch on each clique in enumerated set
//   for (std::vector<NodeID> &clique : curr_cliques) {
//     // add new clique and remove from G
//     reduVCC.addClique(clique);
//     reduVCC.removeVertexSet(clique);
//     // std::cout << "new queue" << std::endl;
//
//     vertex_queue *new_queue = new vertex_queue(G);
//     for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
//
//     // std::cout << "branch" << std::endl;
//     // branch
//     branch_count++;
//     cascading_red_bandr(G, num_fold_cliques, new_queue, partition_config, t);
//
//     // pop branched on clique
//     reduVCC.pop_clique(clique);
//     reduVCC.addVertexSet(clique);
//
//   }
//   // undo number of reductions from reduce
//   R.undoReductions(G, reduVCC); reducer_stack.pop_back();
// }
//
//   // void branch_and_reduce::generate_mis_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t) {
//   //
//   //     if (t.elapsed() > partition_config.solver_time_limit) return;
//   //
//   //     // perform exhaustive reductions
//   //     reducer R(G);
//   //     R.exhaustive_reductions(G, reduVCC); reducer_stack.push_back(R);
//   //     // keep track of total folded cliques to determine current clique cover size
//   //     num_fold_cliques += R.num_fold_cliques;
//   //
//   //     // current size of parital clique cover
//   //     unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;
//   //
//   //     // check exit condition -- kernel is empty
//   //     if (reduVCC.remaining_nodes == 0) {
//   //       // check if we have a better solution
//   //       if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
//   //         // build current parital cover
//   //         // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   //         reduVCC.build_cover(G);
//   //
//   //         // unwind reductions to get full cover
//   //         for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(G, reduVCC);
//   //       }
//   //
//   //       // undo branch's reductions and return
//   //       R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//   //       return;
//   //     }
//   //
//   //     // geneate MIS of kernel using ILS
//   //     graph_access G_p;
//   //     graph_io::readGraphKernel(G_p, reduVCC);
//   //     MISConfig config;
//   //     config.console_log = true;
//   //     config.time_limit = 60;
//   //     config.force_cand = 4;
//   //     ils new_ils;
//   //     new_ils.perform_ils(config, G_p, 1000);
//   //
//   //
//   //     // estimate cover size using partial cover size and MIS of kernel
//   //     unsigned int estimated_cover_size = curr_cover_size + new_ils.solution_size;
//   //     // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
//   //     // prune branch if estimated cover is larger than current best
//   //     if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
//   //       // std::cout << "prune" << std::endl;
//   //       R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//   //       return;
//   //     }
//   //
//   //
//   //     // get next node in kernel with minimum degree
//   //     NodeID next_node = min_deg_node(reduVCC);
//   //
//   //     // enumerate all maximal cliques of next_node sorted by size and MIS
//   //     // std::cout << "enumerate" << std::endl;
//   //     std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node);
//   //     // std::cout << "complete enumerate" << std::endl;
//   //     // branch on each clique in enumerated set
//   //     for (std::vector<NodeID> &clique : curr_cliques) {
//   //       // add new clique and remove from G
//   //       reduVCC.addClique(clique);
//   //       reduVCC.removeVertexSet(clique);
//   //       // std::cout << "branch" << std::endl;
//   //
//   //       // branch
//   //       branch_count++;
//   //       generate_mis_bandr(G, num_fold_cliques, partition_config, t);
//   //
//   //       // pop branched on clique
//   //       reduVCC.pop_clique(clique);
//   //       reduVCC.addVertexSet(clique);
//   //
//   //     }
//   //     // undo number of reductions from reduce
//   //     R.undoReductions(G, reduVCC); reducer_stack.pop_back();
//   //   }
