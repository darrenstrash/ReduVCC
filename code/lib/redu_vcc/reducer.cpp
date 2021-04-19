
#include "reducer.h"


reducer::reducer(graph_access &G) {
  num_reductions = 0;
  num_attempts = 0;
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

void reducer::bruteISO(graph_access &G, redu_vcc &reduVCC, std::vector<unsigned int> &iso_degree) {

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

              iso_degree[pReduction->deg]++;
         }
         num_attempts++;
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
         if (reduVCC.adj_size(v) == 2) num_attempts++;
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
         if (reduVCC.adj_size(v) == 3) num_attempts++;
      } endfor

 }
}

void reducer::bruteDOM(graph_access &G, redu_vcc &reduVCC, std::vector<unsigned int> &dom_degree) {

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

              dom_degree[pReduction->deg]++;

         }
         num_attempts++;
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
  else { delete pReduction; }
  num_attempts++;

}

void reducer::exhaustive_reductions(graph_access &G, redu_vcc &reduVCC,
                                    std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &dom_degree){

  bool new_reduced = true;
  unsigned int curr_reductions = 0;

  while (new_reduced) {
    // std::cout << curr_reductions << std::endl;
    new_reduced = false;
    bruteISO(G, reduVCC, iso_degree);
    bruteD2(G, reduVCC);
    bruteTWIN(G, reduVCC);
    bruteDOM(G, reduVCC, dom_degree);
    bruteCROWN(G, reduVCC);

    if (num_reductions > curr_reductions) {
      curr_reductions = num_reductions;
      new_reduced = true;
    }
  }

}

void reducer::cascading_reductions(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                                   std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &dom_degree){

  while(!queue->empty()) {
    while(!queue->empty()) {
      NodeID v = queue->pop();
      if (!reduVCC.node_status[v]) continue;
      // std::cout << queue->size() << std::endl;

      NodeID u;

      reduction *pReduction = nullptr;

      unsigned int adj_size = reduVCC.adj_size(v);

      if (iso_reduction::validISO(reduVCC, v)) {
        iso_degree[adj_size]++;
        num_attempts++;
        pReduction = new iso_reduction();
      }
      else if (d2_reduction::validD2(reduVCC, v)){
        // if (adj_size <= 5) num_attempts++;
        num_attempts += 2;
        pReduction = new d2_reduction();
      }
      else if (twin_reduction::validTWIN(reduVCC, v, u)){
        // if (adj_size <= 5) num_attempts++;
        if (adj_size == 2) num_attempts++;
        num_attempts+=2;
        pReduction = new twin_reduction();
      }
      else if (dom_reduction::validDOM(reduVCC, v, u)){
        dom_degree[adj_size]++;
        // if (adj_size <= 5) num_attempts++;
        if (adj_size == 2) num_attempts++;
        if (adj_size == 3) num_attempts++;
        num_attempts+=2;
        pReduction = new dom_reduction();
      }
      else {
        num_attempts +=2;
        // if (adj_size <= 5) num_attempts++;
        if (adj_size == 2) num_attempts++;
        if (adj_size == 3) num_attempts++;
        delete pReduction;
        continue;
      }

      pReduction->reduce(G, reduVCC, queue, v, u);
      num_reductions++;
      num_cliques += pReduction-> num_cliques;
      num_fold_cliques += pReduction->num_folded_cliques;
      reduction_stack.push_back(pReduction);

      // if (pReduction->type.compare("iso")){
      //   iso_degree[pReduction->deg]++;
      // }
      // if (pReduction->type.compare("dom")){
      //   dom_degree[pReduction->deg]++;
      // }
    }

    NodeID v; NodeID u;

    reduction *pReduction = new crown_reduction();
    pReduction->reduce(G, reduVCC, queue, v, u);

    if (pReduction->num_cliques > 0) {
      reduction_stack.push_back(pReduction);

      num_reductions++;
      num_cliques += pReduction-> num_cliques;
      // std::cout << "crown cliques: "<< pReduction->num_cliques << std::endl;;
    }
    else { delete pReduction; }
    num_attempts++;
    }
}
