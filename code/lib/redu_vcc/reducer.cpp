
#include "reducer.h"


reducer::reducer(graph_access &G) {
  reduVCC = redu_vcc(G);

  // for_adjList(reduVCC, 0, u) {
  //   std::cout << u << std::endl;
  // } endfor
}

void reducer::buildCover(graph_access &G) {
  reduVCC.build_cover(G);
}

void reducer::unwindReductions(graph_access &G) {

    for (unsigned int i = reduction_stack.size(); i > 0; i--) {

        reduction_stack[i-1]->unfold(G, reduVCC);

    }

    // std::cout << reduVCC.next_cliqueID << std::endl;
}

void reducer::undoReductions(graph_access &G, unsigned int num) {

  for (unsigned int i = 0; i < num; i++) {
      reduction *pReduction = reduction_stack.back();
      reduction_stack.pop_back();

      pReduction->unreduce(G, reduVCC);

      delete pReduction;
  }
}

std::vector<unsigned int> reducer::bruteISO(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;
  unsigned int num_cliques = 0;
  unsigned int num_folded_cliques = 0;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (iso_reduction::validISO(reduVCC, v)){
              reduced++;
              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new iso_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

              num_cliques += pReduction->num_cliques;
              num_folded_cliques += pReduction->num_folded_cliques;
         }
      } endfor

 }
 return {reduced, num_cliques, num_folded_cliques};
}

std::vector<unsigned int> reducer::bruteD2(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;
  unsigned int num_cliques = 0;
  unsigned int num_folded_cliques = 0;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (d2_reduction::validD2(reduVCC, v)){
              reduced++;
              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new d2_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

              num_cliques += pReduction->num_cliques;
              num_folded_cliques += pReduction->num_folded_cliques;

         }
      } endfor

 }

 return {reduced, num_cliques, num_folded_cliques};
}

std::vector<unsigned int> reducer::bruteTWIN(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;
  unsigned int num_cliques = 0;
  unsigned int num_folded_cliques = 0;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (twin_reduction::validTWIN(reduVCC, v, u)){

              // std::cout << "twin" << std::endl;
              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new twin_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

              reduced++;
              num_cliques += pReduction->num_cliques;
              num_folded_cliques += pReduction->num_folded_cliques;

         }
      } endfor

 }

 return {reduced, num_cliques, num_folded_cliques};
}

std::vector<unsigned int> reducer::bruteDOM(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;
  unsigned int num_cliques = 0;
  unsigned int num_folded_cliques = 0;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          // if ( reduced > 12 ) return reduced;

          // std::cout << "starts" << std::endl;

          if (!reduVCC.node_status[v]) { continue;}

          if (dom_reduction::validDOM(reduVCC, v, u)){

              // std::cout<< "valid dom" << std::endl;
              reduced++;
              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new dom_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

              num_cliques += pReduction->num_cliques;
              num_folded_cliques += pReduction->num_folded_cliques;

         }
         // std::cout << "ends" << std::endl;
      } endfor

      // std::cout << "whle agin" << std::endl;

 }

 // std::cout << "returns" << std::endl;
 return {reduced, num_cliques, num_folded_cliques};
}

std::vector<unsigned int> reducer::bruteCROWN(graph_access &G) {

  unsigned int reduced = 0;
  unsigned int num_cliques = 0;
  unsigned int num_folded_cliques = 0;

  NodeID v; NodeID u;

  unsigned int curr_cliqueID = reduVCC.next_cliqueID;

  reduction *pReduction = nullptr;
  pReduction = new crown_reduction();
  pReduction->reduce(G, reduVCC, v, u);
  // reduction_stack.push_back(pReduction);

  if (pReduction->num_cliques > 0) {
    reduction_stack.push_back(pReduction);
    reduced++;
    num_cliques+= pReduction-> num_cliques;
  }
  else { delete pReduction; }

  return {reduced, num_cliques, num_folded_cliques};
  // unsigned int num_crown = reduVCC.next_cliqueID - curr_cliqueID;
  // if (num_crown > 0) {
  //   reduction_stack.push_back(pReduction);
  //   reduced++;
  //   num_crown += num
  // } else { delete pReduction; }
  // std::cout << num_crown << std::endl;
  // return num_crown;
}

void reducer::solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t) {

  unsigned int num_nodes = reduVCC.remaining_nodes;
  if (num_nodes == 0) { return; }

  if (partition_config.solver_time_limit == 0) return;

  reduVCC.buildKernel(G);
  std::vector<std::vector<int>> &kernel_adj_list = reduVCC.kernel_adj_list;
  unsigned long &num_edges = reduVCC.kernel_edges;

  cli *cli_instance;
  cli_instance = new cli(partition_config.seed, partition_config.mis, partition_config.run_type != "upper_bound");
  cli_instance->start_cli(kernel_adj_list, num_nodes, num_edges, t.elapsed(), partition_config.solver_time_limit);

  if (cli_instance->clique_cover.size() != 0){
    // std:: cout << cli_instance->clique_cover.size() << std::endl;
      if (partition_config.run_type != "upper_bound") reduVCC.addKernelCliques(cli_instance->clique_cover);
      chalupa_upper_bound = cli_instance->clique_cover.size();
      chalupa_mis = (unsigned int) cli_instance->final_indset_size;
      std::cout << chalupa_mis << ", " << reduVCC.clique_cover.size() << std::endl;
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
