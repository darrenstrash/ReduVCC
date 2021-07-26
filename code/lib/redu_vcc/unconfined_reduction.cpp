
#include <algorithm>
#include <iostream>
#include <fstream>

#include "unconfined_reduction.h"

// bool unconfined_reduction::nodeDominates(redu_vcc &reduVCC, NodeID &v, NodeID &a){
//
//   std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
//
//   unsigned int i = 0;
//   unsigned int j = 0;
//
//   while (j < adj_list[a].size()){
//       if (!reduVCC.node_status[adj_list[a][j]] || adj_list[a][j] == v){
//           j++;
//           continue;
//       }
//       if (i == adj_list[v].size()){
//           return false;
//       }
//       if (!reduVCC.node_status[adj_list[v][i]] || adj_list[v][i] == a){
//           i++;
//           continue;
//       }
//       if (adj_list[v][i] > adj_list[a][j]){
//           return false;
//       }
//       if (adj_list[v][i] < adj_list[a][j]){
//           i++;
//       }
//       else {
//           i++;
//           j++;
//       }
//   }
//
//   std::vector<NodeID> adj_2 = reduVCC.curr_adj_list(a);
// std::vector<NodeID> adj_1  = reduVCC.curr_adj_list(v);
// adj_2.push_back(a);
// adj_1.push_back(v);
// // std::sort(adj_list[p].begin(), adj_list[p].end());
// // reduVCC.printVectorSet(adj_1);
// // reduVCC.printVectorSet(adj_2);
//
//   // if (!unconfined_reduction::isSubset(reduVCC, adj_2, adj_1)) std::cout << "error" << std::endl;
//   //
//   // reduVCC.printAdjList(v);
//   // reduVCC.printAdjList(a);
//   return true;
// }


bool unconfined_reduction::validUNCONFINED(redu_vcc &reduVCC, NodeID &v, NodeID &u){
    // checks if v is an isolated vertex

    std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
    std::vector<bool> &node_status = reduVCC.node_status;

    std::vector<NodeID> S;
    S.push_back(v);

    std::vector<bool> S_nodes; // marks vertices in S
    S_nodes.assign(reduVCC.num_nodes, false);
    S_nodes[v] = true;

    std::vector<bool> N_S; // marks vertices in N[S]
    N_S.assign(reduVCC.num_nodes, false);
    N_S[v] = true;
    for (NodeID &u_p : adj_list[v]) N_S[u_p] = true;

    while(true) {

      // find u \subset N(S) st. |N(u) \cap S| = 1 and |N(u) \setminus N[S]| is minimized

      std::vector<bool> checked_u; // prevents u from begin checked twice
      checked_u.assign(reduVCC.num_nodes, false);

      bool u_found = false; // if a u st. |N(u) \cap S| = 1 is found

      // unsigned int u_size = reduVCC.num_nodes; // |N(u) \setminus N[S]| for u with minimized value
      std::vector<NodeID> N_u; // N(u) \setminus N[S] for u with minimized value

      for (NodeID &v_p : S) { // iterates through each vertex in S
        for (NodeID &u_p : adj_list[v_p]) {
          if (!node_status[u_p]) continue;
          if (checked_u[u_p]) continue;
          checked_u[u_p] = true; // marks u_p as checked

          std::vector<NodeID> N_u_p;

          unsigned int crossed_S = 0; // counter for |N(u_p) \cap S|
          unsigned int u_p_size = 0; //counter for |N(u_p)|

          for (NodeID &w : adj_list[u]) {
            if (!node_status[w]) continue;
            if (crossed_S > 2) break; // break if |N(u) \cap S| > 1
            if (S_nodes[w]) crossed_S++; // increases counter if w \in S
            // if (!N_S[w]) u_p_size++; // increases counter if w \notin S
            if (!N_S[w]) N_u_p.push_back(w); // increases counter if w \notin S
          }

          if (crossed_S > 2) continue; // |N(u_p) \cap S| != 1, u_p fails condition

          if (N_u_p.size() < N_u.size() || (N_u.size() == 0 && !u_found)) {
            u = u_p;
            N_u = N_u_p;
            u_found = true;
          }
        }

        if (!u_found) return false; // no neighbor such that |N(u_p) \cap S| = 1
        if (N_u.size() == 0) return true;
        if (N_u.size() > 1) return false;

        // add w /in N_u to S, repeat
        NodeID w = N_u[0];
        S.push_back(w);
        S_nodes[w] = true;
        N_S[w] = true;
        for (NodeID &u_p : adj_list[w]) N_S[u_p] = true;

      }



    }

    // if (adj_list[v].size() <= 2) {return false;}
    // if (reduVCC.adj_size(v) <= 2) return false;
    //
    // for (NodeID a : adj_list[v]) {
    //   if (!reduVCC.node_status[a]) continue;
    //   // if (adj_list[v].size() < adj_list[a].size()) { continue; }
    //   if (reduVCC.adj_size(v) < reduVCC.adj_size(a)) continue;
    //
    //   if (unconfined_reduction::nodeDominates(reduVCC, v, a)) {
    //     u = a;
    //     return true;
    //   }
    // }
    // return false;




}


void unconfined_reduction::reduce(  redu_vcc &reduVCC,
                           NodeID &node_v, NodeID &node_u ){

  type = "unconfined";

  v = node_v;
  u = node_u;

  clique.push_back(v);
  clique.push_back(u);

  reduVCC.addClique(clique);
  reduVCC.removeVertexSet(clique);

  // std::cout << "success" << std::endl;
}

void unconfined_reduction::reduce(  redu_vcc &reduVCC, vertex_queue *queue,
                           NodeID &node_v, NodeID &node_u ){
  reduce( reduVCC, node_v, node_u);
  queue->adjust_queue(reduVCC, v);
  queue->adjust_queue(reduVCC, u);
}

void unconfined_reduction::unfold( redu_vcc &reduVCC) {}

void unconfined_reduction::unreduce( redu_vcc &reduVCC){
  reduVCC.pop_clique(clique);
  reduVCC.addVertexSet(clique);
}
