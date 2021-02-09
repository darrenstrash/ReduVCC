
#include "reducer.h"


reducer::reducer(graph_access &G) {

  reduVCC.build(G);
}

void reducer::unwindReductions(graph_access &G) {

    while (reduction_stack.size() > 0){
        reduction *pReduction = reduction_stack.back();
        reduction_stack.pop_back();

        pReduction->unreduce(G, reduVCC);

        delete pReduction;
    }
}

void reducer::bruteISO(graph_access &G) {

  bool vertexReduced = true;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (iso_reduction::validISO(reduVCC, v)){

              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new iso_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);
         }
      } endfor

 }
}

void reducer::bruteD2(graph_access &G) {

  bool vertexReduced = true;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (d2_reduction::validD2(reduVCC, v)){

              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new d2_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

         }
      } endfor

 }
}

void reducer::bruteTWIN(graph_access &G) {

  bool vertexReduced = true;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (twin_reduction::validTWIN(reduVCC, v, u)){

              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new twin_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

         }
      } endfor

 }
}

void reducer::bruteDOM(graph_access &G) {

  bool vertexReduced = true;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (dom_reduction::validDOM(reduVCC, v, u)){

              std::cout<< "valid dom" << std::endl;

              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new dom_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

         }
      } endfor

 }
}

void reducer::bruteCROWN(graph_access &G) {

  NodeID v; NodeID u;

  reduction *pReduction = nullptr;
  pReduction = new crown_reduction();
  pReduction->reduce(G, reduVCC, v, u);
  reduction_stack.push_back(pReduction);
}

void reducer::solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t) {

  unsigned int num_nodes = reduVCC.remaining_nodes;
  if (num_nodes == 0) { return; }

  reduVCC.buildKernel(G);
  std::vector<std::vector<int>> &kernel_adj_list = reduVCC.kernel_adj_list;
  unsigned long &num_edges = reduVCC.kernel_edges;

  cli *cli_instance;
  cli_instance = new cli(partition_config.seed, partition_config.mis);
  cli_instance->start_cli(kernel_adj_list, num_nodes, num_edges, t.elapsed(), partition_config.solver_time_limit);

  if (cli_instance->clique_cover.size() != 0){
      reduVCC.addKernelCliques(cli_instance->clique_cover);
  }
  else {
      std::cout << "Chalupa's algorithm unable to solve in given time." << std::endl;
  }
  delete(cli_instance);

}

void reducer::analyzeGraph(std::string &filename, graph_access &G, timer &t){

    std::cout << filename << ", ";

    std::cout << G.number_of_nodes() << ", ";
    std::cout << G.number_of_edges() << ", ";

    std::cout << reduVCC.remaining_nodes << ", ";

    std::cout << t.elapsed() << ", ";

    std::cout << reduVCC.clique_cover.size() << std::endl;

    reduVCC.validateCover(G);

}



void reducer::branch( std::vector<std::vector<NodeID>> &clique_cover,
                      std::vector<std::vector<NodeID>> &curr_cover,
                      NodeID curr_node
                    ) {

  // chalupa's code on small graphs to test
  // hpc on "medium" -- 50-100
  // http://networkrepository.com/


  // try implementing reductions with this 


  NodeID next_node = curr_node;

  while (!reduVCC.node_status[next_node]) {
    if (next_node >= reduVCC.node_status.size()) {
      if (curr_cover.size() < clique_cover.size() || clique_cover.size() == 0) {
        clique_cover = curr_cover;
      }
      return;
    }
    next_node++;
  }

  std::vector<std::vector<NodeID>> curr_cliques = enumerate(next_node);

  for (std::vector<NodeID> &clique : curr_cliques) {
    curr_cover.push_back(clique);
    for (NodeID x : clique) { reduVCC.node_status[x] = false;}
    branch(clique_cover, curr_cover, next_node);
    curr_cover.pop_back();
    for (NodeID x : clique) { reduVCC.node_status[x] = true;}
  }
}

void reducer::branch_and_bound() {

  std::vector<std::vector<NodeID>> minimal_cover;
  std::vector<std::vector<NodeID>> current_cover;

  branch(minimal_cover, current_cover, 0);

  for (std::vector<NodeID> &clique : minimal_cover) { reduVCC.printVectorSet(clique); }

}


std::vector<std::vector<NodeID>> reducer::enumerate(NodeID v) {

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &node_status = reduVCC.node_status;

  std::vector<NodeID> N_v {v};
  for (NodeID u : adj_list[v]) {
    if (node_status[u]) { N_v.push_back(u); }
  }
  std::sort(N_v.begin(), N_v.end());

  std::vector<std::vector<NodeID>> enum_cliques;
  std::vector<NodeID> curr_clique = {v};

  enumerator(N_v, 0, curr_clique, enum_cliques);

  return enum_cliques;
}

void reducer::enumerator(
                         std::vector<NodeID> &N_v,
                         unsigned int curr_i,
                         std::vector<NodeID> curr_clique,
                         std::vector<std::vector<NodeID>> &enum_cliques) {

  // changes:
  // test with the other code on small graph
  // report only maximal
  // curr_clique as ref

  // to do:
  // instead of enumerate all cliques -- 1 at a time; function -- as soon as you find a maximal you call
  // impliment pivoting -- avoid redundant recursive calls,

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &node_status = reduVCC.node_status;
  std::vector<bool> &scratch1 = reduVCC.scratch1;

  // reduVCC.printVectorSet(curr_clique);

  enum_cliques.push_back(curr_clique);
  NodeID curr_node = N_v[curr_i];

  for (NodeID x : adj_list[curr_node]) {
    if (node_status[x]) { scratch1[x] = true;}
  }
  std::vector<unsigned int> next_i;

  for (unsigned int i = curr_i + 1; i < N_v.size(); i++) {
    NodeID n = N_v[i];
    if (!node_status[n]) { continue; }
    if (scratch1[n]) { next_i.push_back(i); }
  }

  reduVCC.clearScratch(scratch1); // fix this so it takes the set that is marked

  for (unsigned int i : next_i) { // maximal clique if cannot add another vertex
    NodeID n = N_v[i];
    std::vector<NodeID> new_clique = curr_clique; // changes this so it doesn't duplicate
    new_clique.push_back(n);
    enumerator(N_v, i, new_clique, enum_cliques);
  }
}



// ENUMERATION CODE FOR UNDELETED VERTICIES

// void reducer::enumerate(NodeID v) {
//
//   std::vector<NodeID> N_v = reduVCC.adj_list[v];
//   N_v.push_back(v);
//   std::sort(N_v.begin(), N_v.end());
//   std::vector<std::vector<NodeID>> enum_cliques;
//
//   std::vector<NodeID> curr_clique = {v};
//
//   enumerator(N_v, 0, curr_clique, enum_cliques);
// }
//
// void reducer::enumerator(
//                          std::vector<NodeID> &N_v,
//                          unsigned int curr_i,
//                          std::vector<NodeID> curr_clique,
//                          std::vector<std::vector<NodeID>> &enum_cliques) {
//
//   std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
//   std::vector<bool> &scratch1 = reduVCC.scratch1;
//
//   reduVCC.printVectorSet(curr_clique);
//
//     enum_cliques.push_back(curr_clique);
//     NodeID curr_node = N_v[curr_i];
//
//     for (NodeID x : adj_list[curr_node]) { scratch1[x] = true;}
//     std::vector<unsigned int> next_i;
//
//     for (unsigned int i = curr_i + 1; i < N_v.size(); i++) {
//       NodeID n = N_v[i];
//       if (scratch1[n]) { next_i.push_back(i); }
//     }
//
//     reduVCC.clearScratch(scratch1);
//
//     for (unsigned int i : next_i) {
//       NodeID n = N_v[i];
//       std::vector<NodeID> new_clique = curr_clique;
//       new_clique.push_back(n);
//       enumerator(N_v, i, new_clique, enum_cliques);
//     }
//
//   }
