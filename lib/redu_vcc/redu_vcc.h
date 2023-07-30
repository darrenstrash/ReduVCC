
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


class redu_vcc {
private:


/********************** remove **********************************/
  // subgraph mapping
  std::vector<int> old_to_new_map;
  std::vector<NodeID> new_to_old_map;
  void assignMaps();


  // generate adj_list from G
  void generateAdjList(graph_access &G);
  void generateAdjList(redu_vcc *parent, std::vector<NodeID> &child_nodes);

  // initialize class attributes
  void init();

  // finds disconnected component in current graph
  std::vector<NodeID> find_component( std::vector<bool> &visited_nodes, unsigned int &visit_remaining);

  // checks if clique is in original graph G
  bool cliqueInG(graph_access &G, std::vector<NodeID> &clique);


public:

  // mapping from current graph to parent graph
  std::vector<NodeID> self_to_parent;

  unsigned int num_nodes;
  std::vector<std::vector<NodeID>> adj_list;

  std::vector<bool> node_status;  // marks nodes status in G
  std::vector<bool> fold_node;   // marks nodes removed in a fold
  std::vector<bool> merge_node;   // marks nodes removed in merge (edgeBnR)
  unsigned int remaining_nodes;

  std::vector<std::vector<NodeID>> clique_cover;

  unsigned int next_cliqueID; // next cliqueID while reducing
  unsigned int next_solvecliqueID; // next cliqueID while solving
  std::vector<unsigned int> node_clique; // node to cliqueID mapping
  std::vector<unsigned int> solve_node_clique; // node to cliqueID mapping for solving

  // marks vertices in the independent set
  std::vector<bool> node_mis;
  // independent set size
  unsigned int curr_mis;

  std::vector<std::vector<int>> kernel_adj_list; // subgraph adj_list
  unsigned long kernel_edges; // num subgraph edges

  std::vector<bool> scratch1;
  std::vector<bool> scratch2;


  redu_vcc() {};
  redu_vcc(graph_access &G);
  redu_vcc(redu_vcc *parent, std::vector<NodeID> &child_nodes);
  redu_vcc(graph_access &G, PartitionConfig &partition_config);
  virtual ~redu_vcc() {};

  // graph methods
  std::vector<redu_vcc> decompose();

  void removeVertex(NodeID v);
  void addVertex(NodeID v);
  void removeVertexSet(std::vector<NodeID> &S);
  void addVertexSet(std::vector<NodeID> &S);

  unsigned int adj_size(NodeID v);
  std::vector<NodeID> curr_adj_list(NodeID v);

  void printAdjList();
  void printAdjList(NodeID v);
  void printNeighborhood(NodeID v);
  void printVectorSet(std::vector<NodeID> S);

  // cover methods
  void addCliquesToParent(redu_vcc &parent);

  void addClique(std::vector<NodeID> &clique);
  void addCliqueToCover(std::vector<NodeID> &clique);
  void pop_clique(std::vector<NodeID> &clique);

  unsigned int getCliqueID(NodeID &v);
  std::vector<NodeID> getClique(NodeID &v);
  void replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique);

  void build_cover();
  bool validateCover(graph_access &G);
  void analyzeGraph(std::string &filename, graph_access &G, timer &t, bool const validate_cover = true);

  // mis methods
  void getMIS(std::string file);

  // integer kernel methods
  void buildKernel();
  void solveKernel(PartitionConfig &partition_config, timer &t, double &time_to_solution, std::size_t clique_cover_offset);
  void addKernelCliques(std::vector<std::vector<int>> &clique_set);

  void addCrownCliques(std::vector<std::vector<NodeID>> &crown_cliques, std::vector<std::vector<int>> &clique_set);

};



#endif
