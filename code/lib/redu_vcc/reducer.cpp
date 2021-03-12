
#include "reducer.h"


reducer::reducer(graph_access &G) {
  num_reductions = 0;
  num_cliques = 0;
  num_fold_cliques = 0;
}


void reducer::unwindReductions(graph_access &G, redu_vcc &reduVCC) {

    for (unsigned int i = reduction_stack.size(); i > 0; i--) {

        reduction_stack[i-1]->unfold(G, reduVCC);

    }

    // std::cout << reduVCC.next_cliqueID << std::endl;
}

void reducer::undoReductions(graph_access &G, redu_vcc &reduVCC) {

  while (reduction_stack.size() != 0) {
        reduction *pReduction = reduction_stack.back();
        reduction_stack.pop_back();

        pReduction->unreduce(G, reduVCC);

        delete pReduction;
  }

  // for (unsigned int i = 0; i < num; i++) {
  //     reduction *pReduction = reduction_stack.back();
  //     reduction_stack.pop_back();
  //
  //     pReduction->unreduce(G, reduVCC);
  //
  //     delete pReduction;
  // }
}

void reducer::bruteISO(graph_access &G, redu_vcc &reduVCC) {

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

              num_reductions++;
              num_cliques++;;
         }
      } endfor

 }
}

void reducer::bruteD2(graph_access &G, redu_vcc &reduVCC) {

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

              num_reductions++;
              num_fold_cliques++;

         }
      } endfor

 }
}

void reducer::bruteTWIN(graph_access &G, redu_vcc &reduVCC) {

  bool vertexReduced = true;

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

              num_reductions++;
              num_cliques += pReduction->num_cliques;
              num_fold_cliques += pReduction->num_folded_cliques;

         }
      } endfor

 }
}

void reducer::bruteDOM(graph_access &G, redu_vcc &reduVCC) {

  bool vertexReduced = true;

  while (vertexReduced){
      vertexReduced = false;

      forall_nodes(G, v){
          NodeID u;

          if (!reduVCC.node_status[v]) { continue;}

          if (dom_reduction::validDOM(reduVCC, v, u)){

              // std::cout<< "valid dom" << std::endl;
              vertexReduced = true;
              reduction *pReduction = nullptr;
              pReduction = new dom_reduction();
              pReduction->reduce(G, reduVCC, v, u);
              reduction_stack.push_back(pReduction);

              num_reductions++;

         }
      } endfor

 }

}

void reducer::bruteCROWN(graph_access &G, redu_vcc &reduVCC) {

  NodeID v; NodeID u;

  unsigned int curr_cliqueID = reduVCC.next_cliqueID;

  reduction *pReduction = nullptr;
  pReduction = new crown_reduction();
  pReduction->reduce(G, reduVCC, v, u);
  // reduction_stack.push_back(pReduction);

  if (pReduction->num_cliques > 0) {
    reduction_stack.push_back(pReduction);

    num_reductions++;
    num_cliques += pReduction-> num_cliques;
  }
  else {
    delete pReduction; }


}

void reducer::exhaustive_reductions(graph_access &G, redu_vcc &reduVCC){

  bool new_reduced = true;
  unsigned int curr_reductions = 0;

  while (new_reduced) {
    // std::cout << curr_reductions << std::endl;
    new_reduced = false;
    bruteISO(G, reduVCC);
    // std::cout << "finish iso" << std::endl;
    bruteD2(G, reduVCC);
    // std::cout << "finish d2" << std::endl;
    bruteTWIN(G, reduVCC);
    // std::cout << "finish twin" << std::endl;
    bruteDOM(G, reduVCC);
    // std::cout << "finish dom" << std::endl;
    bruteCROWN(G, reduVCC);
    // std::cout << "finish crown" << std::endl;


    if (num_reductions > curr_reductions) {
      curr_reductions = num_reductions;
      new_reduced = true;
    }
  }

}
