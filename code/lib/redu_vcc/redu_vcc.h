
/******************************************************************************
 * redu_vcc.h
 * *
 *
 *****************************************************************************/

#ifndef REDU_VCC
#define REDU_VCC

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_structure.h"


class redu_vcc : public redu_structure {

private:

  std::vector<NodeID> parent_to_self_map;
  std::vector<NodeID> self_to_parent_map;

  void generateAdjList(graph_access &G);
  void generateAdjList(redu_vcc* parent);
  void subgraph_map(std::vector<NodeID> &subgraph_nodes);

  // void generateReducedAdjList(std::vector<int> &node_to_int_map, std::vector<NodeID> &int_to_node_map,
  //                             std::vector<std::vector<int>> &int_adj_list, unsigned long &edges_count);
  // void addSolveCliques(std::vector<NodeID> &int_to_node_map, std::vector<std::vector<int>> &solve_cliques);
  //


public:

  redu_vcc() {};
  redu_vcc(graph_access &G, PartitionConfig &partition_config);
  redu_vcc(redu_vcc *parent, std::vector<NodeID> &subgraph_nodes);
  virtual ~redu_vcc() {};

  void getMIS(std::string file);

  std::vector<redu_vcc> decompose_components();
  std::vector<NodeID> find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining);


  void merge_covers(redu_vcc &parent);

  void analyzeGraph(std::string &filename, graph_access &G, timer &t);
  void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);


};

#endif
