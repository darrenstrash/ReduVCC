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
  // prune_type = "KaMIS";
  // if (partition_config.run_type == "KaMIS") return;
  redu_type = "cascading";
  if (partition_config.run_type == "cascading") return;
}

branch_and_reduce::branch_and_reduce(graph_access &G, PartitionConfig &partition_config) {

  construct_run(partition_config);

  if (prune_type == "ReduMIS") reduVCC = redu_vcc(G, partition_config);
  else reduVCC = redu_vcc(G);

  branch_count = 0;
  iso_degree.assign(G.number_of_nodes(), 0);
  dom_degree.assign(G.number_of_nodes(), 0);
  num_reductions = 0;
  num_attempts = 0;

  // config.time_limit = 60;
  // config.force_cand = 4;
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

// std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x, std::vector<bool> &indset) {
std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(NodeID x) {

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(x);
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

void branch_and_reduce::reduce(graph_access &G, reducer &R, unsigned int &num_fold_cliques, vertex_queue *queue) {

    if (queue == nullptr || queue->empty()) {
      R.exhaustive_reductions(G, reduVCC, iso_degree, dom_degree);
    }
    else {
      R.cascading_reductions(G, reduVCC, queue, iso_degree, dom_degree);
    }
    reducer_stack.push_back(R);

    num_reductions += R.num_reductions;
    num_attempts += R.num_attempts;

    // keep track of total folded cliques to determine current clique cover size
    num_fold_cliques += R.num_fold_cliques;

}

bool branch_and_reduce::prune(unsigned int &curr_cover_size) {

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
    else { // prune_type == "ReduMIS"
      estimated_cover_size = curr_cover_size + reduVCC.curr_mis;
    }

    // std::cout << "est cover: " << estimated_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    // prune branch if estimated cover is larger than current best
    if (reduVCC.clique_cover.size() != 0 && estimated_cover_size >= reduVCC.clique_cover.size()) {
      // std::cout << "prune" << std::endl;
      return true;
    }

    return false;

}

NodeID branch_and_reduce::min_deg_node() {

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

std::vector<std::vector<NodeID>> branch_and_reduce::enum_vertex(NodeID &v) {

  if (enum_type == "sort_enum") return sorted_enumerate(v);
  else return enumerate(v);
}

vertex_queue* branch_and_reduce::construct_queue(graph_access &G, std::vector<NodeID> &clique) {

  vertex_queue *new_queue = nullptr;
  if (redu_type == "exhaustive") return new_queue;

  new_queue = new vertex_queue(G);
  for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);
  return new_queue;

}

NodeID branch_and_reduce::nextNode(){

  if (next_node_type == "small_deg") {
    return min_deg_node();
  }
  NodeID next_node = 0;
  while (!reduVCC.node_status[next_node]) next_node++;
  return next_node;

}

void branch_and_reduce::bandr( graph_access &G, unsigned int num_fold_cliques,
                               vertex_queue *queue,
                               PartitionConfig &partition_config, timer &t) {

  if (t.elapsed() > partition_config.solver_time_limit) return;

  reducer R(G, partition_config.iso_limit);
  reduce(G, R, num_fold_cliques, queue);
  delete queue;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;

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


  if (prune(curr_cover_size)) {
    R.undoReductions(G, reduVCC); reducer_stack.pop_back();
    return;
  }

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
  NodeID next_node = nextNode();

  // enumerate all maximal cliques of next_node sorted by size and MIS
  // std::cout << "enumerate" << std::endl;
  // std::vector<std::vector<NodeID>> curr_cliques = sorted_enumerate(next_node, reduVCC.node_mis);
  std::vector<std::vector<NodeID>> curr_cliques = enum_vertex(next_node);

  // std::cout << "complete enumerate" << std::endl;
  // branch on each clique in enumerated set
  for (std::vector<NodeID> &clique : curr_cliques) {
    // add new clique and remove from G
    reduVCC.addClique(clique);
    reduVCC.removeVertexSet(clique);
    // std::cout << "new queue" << std::endl;

    vertex_queue *new_queue = construct_queue(G, clique);
    // vertex_queue *new_queue = new vertex_queue(G);
    // for (NodeID a : clique) new_queue->adjust_queue(reduVCC, a);

    // std::cout << "branch" << std::endl;
    // branch
    branch_count++;
    bandr(G, num_fold_cliques, new_queue, partition_config, t);

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(G, reduVCC); reducer_stack.pop_back();
}

std::vector<NodeID> branch_and_reduce::find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining) {

  std::vector<NodeID> current_nodes;
  std::vector<NodeID> queue;

  NodeID v = 0;
  while (visited_nodes[v]) v ++;
  visited_nodes[v] = true; visit_remaining--;
  current_nodes.push_back(v);

  queue.push_back(v);

  while (!queue.size() == 0) {

    v = queue.front();
    queue.erase(queue.begin());


    for (NodeID u : reduVCC.adj_list[v]) {
      if (!visited_nodes[u]) {
        visited_nodes[u] = true; visit_remaining--;
        current_nodes.push_back(u);
        queue.push_back(u);
      }
    }

  }

  for (NodeID a : current_nodes) {
    std::cout << a << ", ";
  }
  std::cout << std::endl;

  return current_nodes;
}

void branch_and_reduce::unconnected_components() {

  std::vector<bool> visited_nodes;
  for (bool n : reduVCC.node_status) visited_nodes.push_back(!n);
  unsigned int visit_remaining = reduVCC.node_status.size();

  std::vector<redu_vcc> components;

  std::vector<NodeID> component_nodes = find_component(visited_nodes, visit_remaining);
  if (visit_remaining == 0) {
    components.push_back(reduVCC)
    return components
  }

  redu_vcc new_component(reduVCC, component_nodes);
  component_nodes.push_back(new_component);

  while (visit_remaining > 0) {
    component_nodes = find_component(visited_nodes, visit_remaining);
    new_component = redu_vcc(reduVCC, component_nodes)
    // redu_vcc new_component (reduVCC, component_nodes)
    componens.push_back(new_component)
  }

  return components;

}
