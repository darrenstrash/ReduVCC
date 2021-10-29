#include <iostream>
#include <fstream>

#include "b_and_r.h"


void branch_and_reduce::construct_run(PartitionConfig &partition_config) {

  redu_type = partition_config.redu_type;
  // redu_type = partition_config.run_type;
  prune_type = partition_config.prune_type;
  if (partition_config.run_type == "edge_bnr") return;

  next_node_type = "small_deg";
  enum_type = "sorted_enum";

}

branch_and_reduce::branch_and_reduce(graph_access &G, redu_vcc &reduVCC, PartitionConfig &partition_config) {

  construct_run(partition_config);

  if (prune_type == "reduMIS") reduVCC = redu_vcc(G, partition_config);
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
    }
  }

  for (NodeID y : N_pivot) reduVCC.scratch1[y] = false;


  std::vector<NodeID> pivot_consider {};
  for (NodeID x : consider_nodes) {
    if (reduVCC.scratch1[x]) pivot_consider.push_back(x);
  }

  for (NodeID y : consider_nodes) reduVCC.scratch1[y] = false;

  while (!pivot_consider.empty()) {

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

  std::vector<std::vector<NodeID>> branch_and_reduce::sorted_enumerate(redu_vcc &reduVCC, NodeID x) {

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(reduVCC, x);

  // sort enumerated cliques
  std::vector<unsigned int> curr_cliques_indices;
  std::vector<unsigned int> curr_clique_sizes;

  for (int i =0; i < curr_cliques.size(); i++) {

    curr_cliques_indices.push_back(i);
    curr_clique_sizes.push_back(curr_cliques[i].size());

  }

  //order cliques
  std::sort(curr_cliques_indices.begin(), curr_cliques_indices.end(),
    [curr_clique_sizes](unsigned int i, unsigned int j) {

          return curr_clique_sizes[i] > curr_clique_sizes[j];
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
    // else if (prune_type == "ils") {
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
    else if (prune_type == "sigmod_linear") {
      // geneate MIS of kernel using Sigmod MIS
      Graph mis_G;
      // std::cout << "begin graph trans" << std::endl;
      mis_G.read_graph(reduVCC);
      // std::cout << "graph read"<< reduVCC.remaining_nodes << std::endl;
      // reduVCC.printAdjList();
      unsigned int res_mis = mis_G.degree_two_kernal_and_remove_max_degree_without_contraction();
      // std::cout << "mis computed" << std::endl;
      estimated_cover_size = curr_cover_size + res_mis;
    }
    else if (prune_type == "sigmod_nearlinear") {
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
    if (t.elapsed() > partition_config.solver_time_limit) return false;

    // solve child
    // std::cout << "next child" << std::endl;
    branch_and_reduce B_child(child, partition_config);
    vertex_queue *queue = nullptr;
    if (partition_config.run_type == "cascading") queue = new vertex_queue(child);
    // child.printAdjList();
    if (!B_child.bandr(child, 0, queue, partition_config, t))
        return false;
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

  if (t.elapsed() > partition_config.solver_time_limit) return false;
  // std::cout << "cover size: " << cover_size << std::endl;

  // check clique cover size
  if (reduVCC.clique_cover.size() == 0 || cover_size < reduVCC.clique_cover.size()) {
    // for each child, add cliques to cover


    //std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    reduVCC.build_cover();
    // std::cout << "builds?" << std::endl;
    for (redu_vcc &child : children) child.addCliquesToParent(reduVCC);

    // unwind reductions to get full cover
    for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
  }

  return true;
}

bool branch_and_reduce::bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                               vertex_queue *queue,
                               PartitionConfig &partition_config, timer &t) {

  if (t.elapsed() > partition_config.solver_time_limit) return false;

  reducer R(partition_config.iso_limit);
  reduce(reduVCC, R, num_fold_cliques, queue);
  delete queue;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      reduVCC.build_cover();

      // unwind reductions to get full cover
      for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }


  if (prune(reduVCC, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }

  if (reduVCC.remaining_nodes > partition_config.decompose_limit && decompose(reduVCC, partition_config, t, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }


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
    if (!bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t))
        return false;

    // pop branched on clique
    reduVCC.pop_clique(clique);
    reduVCC.addVertexSet(clique);

  }
  // undo number of reductions from reduce
  R.undoReductions(reduVCC); reducer_stack.pop_back();
  return true;
}

bool branch_and_reduce::edge_decompose(redu_vcc &reduVCC, PartitionConfig &partition_config, timer &t,
                                  unsigned int curr_cover_size) {

  // std::cout << "pre decompose remaining nodes: " << reduVCC.remaining_nodes << std::endl;
  // std::cout << "begin decompose" << std::endl;
  std::vector<redu_vcc> children = reduVCC.decompose();

  // std::cout << "post decompose remaining nodes: " << reduVCC.remaining_nodes << std::endl;
  // std::cout << "end decopose" << std::endl;

  if (children.empty()) return false;
  // std::cout << "num_children: "<< children.size() << std::endl;

  unsigned int cover_size = curr_cover_size;
  // std::cout << "curr cover: " << cover_size << std::endl;

  decompose_count += children.size();

  for (redu_vcc &child : children) {
    if (t.elapsed() > partition_config.solver_time_limit) return false;

    // solve child
    // std::cout << "next child " << child.num_nodes << std::endl;
    branch_and_reduce B_child(child, partition_config);
    vertex_queue *queue = nullptr;
    if (partition_config.run_type == "cascading") queue = new vertex_queue(child);
    // std::cout << "child map: [";
    // for (unsigned int i = 0; i < child.self_to_parent.size(); i++) {
    //   std::cout << "(" << i << ", " << child.new_to_old_map[i] << "), ";
    // }
    // std::cout << "]" << std::endl;
    // if (partition_config.run_type == "cascading") queue = new vertex_queue(child);
    // child.printAdjList();
    // std::cout << "new child branch" << std::endl;
    if (!B_child.edge_bandr(child, 0, queue, partition_config, t, 0))
        return false;
    // if (child.num_nodes != child.remaining_nodes) std::cout << "unmatched nodes sizes" << std::endl;
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

  // std::cout << "post solve children remaining nodes: " << reduVCC.remaining_nodes << std::endl;
  //
  // std::cout << "rturn children" << std::endl;

  if (t.elapsed() > partition_config.solver_time_limit) return false;
  // std::cout << "cover size: " << cover_size << std::endl;

  // check clique cover size
  if (reduVCC.clique_cover.size() == 0 || cover_size < reduVCC.clique_cover.size()) {
    // for each child, add cliques to cover

    // std::cout << "merge children smaller cover: " << cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
    reduVCC.build_cover();
    // std::cout << "builds?" << std::endl;
    for (redu_vcc &child : children) child.addCliquesToParent(reduVCC);

    // unwind reductions to get full cover
    // for (unsigned int i = reducer_stack.size(); i > 0; i--) reducer_stack[i-1].unwindReductions(reduVCC);

    for (unsigned int i = reducer_stack.size(); i > 0; i--) {
      reducer_stack[i-1].unwindReductions(reduVCC);

      if (i > 1) {
        std::vector<NodeID> &merged_edge = edge_stack[i-2];
        if (merged_edge.size() == 0) continue;

        NodeID &e_v = merged_edge[0];
        NodeID &e_u = merged_edge[1];

        unsigned int solveCliqueID = reduVCC.solve_node_clique[e_v];
        reduVCC.solve_node_clique[e_u] = solveCliqueID;
        reduVCC.clique_cover[solveCliqueID].push_back(e_u);

        std::sort(reduVCC.clique_cover[solveCliqueID].begin(), reduVCC.clique_cover[solveCliqueID].end());
      }
    }
  }

  return true;
}

void branch_and_reduce::buildcover_edge ( redu_vcc &reduVCC) {

  reduVCC.build_cover();

  // unwind reductions to get full cover
  for (unsigned int i = reducer_stack.size(); i > 0; i--) {
    reducer_stack[i-1].unwindReductions(reduVCC);

    if (i > 1) {
      std::vector<NodeID> &merged_edge = edge_stack[i-2];
      if (merged_edge.size() == 0) continue;

      NodeID &e_v = merged_edge[0];
      NodeID &e_u = merged_edge[1];

      unsigned int solveCliqueID = reduVCC.solve_node_clique[e_v];

      reduVCC.solve_node_clique[e_u] = solveCliqueID;
      reduVCC.clique_cover[solveCliqueID].push_back(e_u);

      std::sort(reduVCC.clique_cover[solveCliqueID].begin(), reduVCC.clique_cover[solveCliqueID].end());
    }
  }
}


bool branch_and_reduce::check_adj(redu_vcc &reduVCC) {

  for (NodeID v = 0; v < reduVCC.num_nodes; v++) {
    if (!reduVCC.node_status[v]) continue;
    std::cout << v << " " << std::endl;
    for (NodeID u : reduVCC.adj_list[v]) {
      if (!reduVCC.node_status[u]) continue;
      std::cout << u << " " << std::endl;
      bool edge_found = false;
      for (NodeID w : reduVCC.adj_list[u]) {
        std::cout << w << " " << std::endl;
        if (v == w) edge_found = true;
      }
      if (! edge_found) return false;
    }
  }
  return true;
}

bool branch_and_reduce::edge_bandr( redu_vcc &reduVCC, unsigned int num_fold_cliques,
                                    vertex_queue *queue, PartitionConfig &partition_config, timer &t,
                                    NodeID curr_node ) {
  /* Branches on edges rather than vertices -- an edge is in a clique cover or not */

  if (t.elapsed() > partition_config.solver_time_limit) return false;

  reducer R(partition_config.iso_limit);
  reduce(reduVCC, R, num_fold_cliques, queue);
  delete queue;

  // current size of parital clique cover
  unsigned int curr_cover_size = reduVCC.next_cliqueID + num_fold_cliques;

  // check exit condition -- kernel is empty
  if (reduVCC.remaining_nodes == 0) {
    // check if we have a better solution
    if (reduVCC.clique_cover.size() == 0 || curr_cover_size < reduVCC.clique_cover.size()) {
      // build current parital cover
      // std::cout << "smaller cover: " << curr_cover_size << ", " << reduVCC.clique_cover.size() << std::endl;
      buildcover_edge(reduVCC);
    }

    // undo branch's reductions and return
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }

  if (prune(reduVCC, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }

  // If there disconnected components, solve individual components
  if (reduVCC.remaining_nodes > partition_config.decompose_limit && edge_decompose(reduVCC, partition_config, t, curr_cover_size)) {
    R.undoReductions(reduVCC); reducer_stack.pop_back();
    return true;
  }


  // get next edge
  curr_node = 0;
  NodeID edge_node;
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

unsigned int branch_num = branch_count;

// edge in some clique
std::vector<NodeID> edge {curr_node, edge_node};

std::vector<NodeID> N_curr;
for (NodeID v : reduVCC.adj_list[curr_node]) N_curr.push_back(v); // save original neighborhood
reduVCC.adj_list[curr_node].clear();

for (NodeID v : reduVCC.adj_list[edge_node])  reduVCC.scratch1[v] = true;
for (NodeID v : N_curr) {
  if (reduVCC.scratch1[v]) { // in both add to common neighborhood
    reduVCC.adj_list[curr_node].push_back(v);
  }
  else { // otherwise remove edge -- must remove curr_node from adj of v
    for (unsigned int i = 0; i < reduVCC.adj_list[v].size(); i++) {
      if (reduVCC.adj_list[v][i] == curr_node) {
        reduVCC.adj_list[v].erase(reduVCC.adj_list[v].begin() + i);
        break;
      }
    }
  }
}
for (NodeID v : reduVCC.adj_list[edge_node])  reduVCC.scratch1[v] = false; reduVCC.scratch1[curr_node] = false;

reduVCC.removeVertex(edge_node);
reduVCC.merge_node[edge_node] = true;
edge_stack.push_back(edge);

// branch on edge
vertex_queue *new_queue = nullptr;
if (redu_type != "exhaustive") {
  new_queue = new vertex_queue(reduVCC);
  new_queue->push(curr_node);
  for (NodeID a : N_curr) if (reduVCC.node_status[a]) new_queue->push(a);
  new_queue->adjust_queue(reduVCC, edge_node);
}
branch_count++;
if (!edge_bandr(reduVCC, num_fold_cliques, new_queue, partition_config, t, curr_node)) return false;

for (NodeID v : reduVCC.adj_list[curr_node])  reduVCC.scratch1[v] = true;
for (NodeID v : N_curr) {
  if (!reduVCC.scratch1[v]) {
    reduVCC.adj_list[v].push_back(curr_node);
    std::sort(reduVCC.adj_list[v].begin(), reduVCC.adj_list[v].end());
    reduVCC.adj_list[curr_node].push_back(v);
  }
}
std::sort(reduVCC.adj_list[curr_node].begin(), reduVCC.adj_list[curr_node].end());
for (NodeID v : reduVCC.adj_list[curr_node])  reduVCC.scratch1[v] = false;

reduVCC.addVertex(edge_node); // add edge node back into graph
reduVCC.merge_node[edge_node] = false;

// edge not in some clique

edge_stack[edge_stack.size()-1].clear();

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

vertex_queue *new_queue2 = nullptr;
if (redu_type != "exhaustive") {
  new_queue2 = new vertex_queue(reduVCC);
  new_queue2->push(curr_node);
  new_queue2->push(edge_node);
}
// branch on no edge
branch_count++;
if (!edge_bandr(reduVCC, num_fold_cliques, new_queue2, partition_config, t, curr_node)) return false;

reduVCC.adj_list[curr_node].push_back(edge_node);
std::sort(reduVCC.adj_list[curr_node].begin(), reduVCC.adj_list[curr_node].end());
reduVCC.adj_list[edge_node].push_back(curr_node);
std::sort(reduVCC.adj_list[edge_node].begin(), reduVCC.adj_list[edge_node].end());

edge_stack.pop_back();

  R.undoReductions(reduVCC); reducer_stack.pop_back();
  return true;
}
