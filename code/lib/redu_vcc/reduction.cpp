
#include <iostream>
#include <fstream>

#include "reduction.h"

bool reduction::isSubset(redu_vcc &reduVCC, std::vector<NodeID> &A, std::vector<NodeID> &B){
  // tests if A is a subset of B

  std::vector<bool> &scratch1 = reduVCC.scratch1;

  for (NodeID v : B) {scratch1[v] = true;}
  for (NodeID v : A) {
      // if (!reduVCC.node_status[v]) continue;
      if (!scratch1[v]) {
          for (NodeID v : B) {scratch1[v] = false;}
          return false;
      }
  }
  for (NodeID v : B) {scratch1[v] = false;}
  return true;
}

void reduction::merge_neighborhoods(redu_vcc &reduVCC, std::vector<NodeID> &disjoint,
                                    std::vector<NodeID> &N_b, NodeID &a, NodeID &b) {
  // merge N[b] into N[a]
  // also marks nodes merged into N_a, not in original N_a
  // constructs N_b

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
  // ensures that there are no edges between N[a] \setminus N[b] and N[b] \setminus N[a]

  std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
  std::vector<bool> &scratch1 = reduVCC.scratch1;
  std::vector<bool> &scratch2 = reduVCC.scratch2;

  // scratch1 = N[a] / N[b]
  // scratch 2 = N[b] / N[a]
  for (NodeID x : adj_list[a]) { scratch1[x] = true;}
  for (NodeID x : adj_list[b]) {scratch2[x] = true; scratch1[x] = false;}
  for (NodeID x : adj_list[a]) {scratch2[x] = false;}

  for (NodeID x : adj_list[a]){
      // if (!reduVCC.node_status[x]) { continue; }
      if (!reduVCC.node_status[x]) break;

      if (!scratch1[x]) {continue;}
      for (NodeID y : adj_list[x]) {
          // if (!reduVCC.node_status[y]) { continue; }
          if (!reduVCC.node_status[y]) break;
          if (scratch2[y]) {
            // reduVCC.clearScratch(scratch1);
            // reduVCC.clearScratch(scratch2);
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
