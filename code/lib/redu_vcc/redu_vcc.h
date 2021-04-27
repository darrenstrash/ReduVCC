
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
#include "reducer.h"


class redu_vcc : public redu_structure {

private:

  std::vector<reducer> reducer_stack;

  std::vector<redu_vcc> children;

  std::vector<NodeID> self_to_parent_map;
  std::vector<NodeID> parent_to_self_map;

  std::vector<std::vector<NodeID>> generateAdjList(graph_access &G);
  unsigned int subgraph_map(std::vector<NodeID> &subgraph_nodes);
  std::vector<std::vector<NodeID>> generateAdjList(redu_vcc *parent, unsigned int &num_nodes);

  std::vector<NodeID> find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining);

public:
  // std::vector<redu_vcc> children;

  redu_vcc() {};
  redu_vcc(graph_access &G);
  redu_vcc(redu_vcc *parent, std::vector<NodeID> &subgraph_nodes);
  virtual ~redu_vcc() {};

  void decompose_components();


    // redu_vcc() {};
    // redu_vcc(graph_access &G) : redu_structure(G) {};
    // redu_vcc(graph_access &G, PartitionConfig &partition_config);
    // virtual ~redu_vcc() {};
    //
    // void getMIS(std::string file);
    //
    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
    // void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);


};

#endif
