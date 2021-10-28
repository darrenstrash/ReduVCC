
#include <iostream>
#include <fstream>

#include "reduction.h"

bool reduction::isSubset(redu_vcc &reduVCC, std::vector<NodeID> &A, std::vector<NodeID> &B){
  /* Tests if A is a subset of B -- note does not handle if elements in A or B are not in G */

  std::vector<bool> &scratch1 = reduVCC.scratch1;

  for (NodeID v : B) scratch1[v] = true;

  for (NodeID v : A) {
      if (!scratch1[v]) {
          for (NodeID v : B) scratch1[v] = false;
          return false;
      }
  }

  for (NodeID v : B) scratch1[v] = false;

  return true;
}

void reduction::merge_neighborhoods(redu_vcc &reduVCC, std::vector<NodeID> &disjoint,
                                    std::vector<NodeID> &N_b, NodeID &a, NodeID &b) {
  /* Merges N[b] into N[a], saves current N[b] (removes vertices not in G),
     and saves N(b) / N(a) */

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;

  for (NodeID p : adj_list[a]) reduVCC.scratch1[p] = true;

  for (NodeID p : adj_list[b]){
      if (!reduVCC.node_status[p]) continue;
      N_b.push_back(p);

      if (p == a) continue;
      if (reduVCC.scratch1[p]) continue;

      adj_list[a].push_back(p);
      adj_list[p].push_back(a);
      disjoint.push_back(p);

      std::sort(adj_list[p].begin(), adj_list[p].end());
  }
  std::sort(adj_list[a].begin(), adj_list[a].end());

  for (NodeID p : adj_list[a]) reduVCC.scratch1[p] = false;
}

bool reduction::uncrossedSets(redu_vcc &reduVCC, NodeID &a, NodeID &b) {
  /* Determines if a and b are crossing independent -- there are no edges
     between N(a) / N(b) and N(b) / N(a) */

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;
  std::vector<bool> &scratch2 = reduVCC.scratch2;

  // scratch1 = N[a] / N[b]
  // scratch 2 = N[b] / N[a]
  for (NodeID x : adj_list[a]) scratch1[x] = true;
  for (NodeID x : adj_list[b]) {scratch2[x] = true; scratch1[x] = false;}
  for (NodeID x : adj_list[a]) scratch2[x] = false;

  for (NodeID x : adj_list[a]){
      if (!reduVCC.node_status[x]) continue;
      if (!scratch1[x]) continue;

      for (NodeID y : adj_list[x]) {
          if (!reduVCC.node_status[y]) continue;
          if (scratch2[y]) {
            for (NodeID x : adj_list[a]) scratch1[x] = false;
            for (NodeID x : adj_list[b]) scratch2[x] = false;
            return false;
          }
      }
  }

  for (NodeID x : adj_list[a]) scratch1[x] = false;
  for (NodeID x : adj_list[b]) scratch2[x] = false;
  return true;

}
