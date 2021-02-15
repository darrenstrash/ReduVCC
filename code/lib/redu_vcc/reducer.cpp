
#include "reducer.h"


void reducer::init(graph_access &G) {

  reduVCC.build(G);
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

unsigned int reducer::bruteISO(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;

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
         }
      } endfor

 }
 return reduced;
}

unsigned int reducer::bruteD2(graph_access &G) {

  bool vertexReduced = true;
  unsigned int reduced = 0;

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

         }
      } endfor

 }

 return reduced;
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
