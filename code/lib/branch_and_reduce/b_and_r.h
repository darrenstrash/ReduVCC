
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef B_AND_R
#define B_AND_R

#include <algorithm>
#include <functional>

#include "data_structure/graph_access.h"
#include "redu_vcc/reducer.h"
#include "redu_vcc/redu_vcc.h"

class branch_and_reduce {
private:
  std::vector<bool> visited_nodes;

  unsigned int mis;
  std::vector<bool> node_mis;

  std::vector<reducer> reducer_stack;
  unsigned int num_reductions;

  void pivot_enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                        std::vector<NodeID> &consider_nodes,
                        std::vector<NodeID> &curr_clique,
                        std::vector<NodeID> &excluded_nodes);


public:

  redu_vcc reduVCC;

  unsigned int branch_count;

  branch_and_reduce(graph_access &G, PartitionConfig &partition_config);
  virtual ~branch_and_reduce() {};

  std::vector<std::vector<NodeID>> enumerate(NodeID v);
  std::vector<std::vector<NodeID>> sorted_enumerate(NodeID x);

  void brute_bandr(graph_access &G, PartitionConfig &partition_config, timer &t, unsigned int num_folded_cliques);
  void mis_bound_bandr(graph_access &G, PartitionConfig &partition_config, timer &t, unsigned int num_folded_cliques);
  void small_degree_bandr(graph_access &G, PartitionConfig &partition_config, timer &t, unsigned int num_folded_cliques);
  void sorted_enum_bandr(graph_access &G, PartitionConfig &partition_config, timer &t, unsigned int num_folded_cliques);

  void analyzeGraph(std::string &filename, graph_access &G, timer &t) {reduVCC.analyzeGraph(filename, G, t);};
};
#endif
