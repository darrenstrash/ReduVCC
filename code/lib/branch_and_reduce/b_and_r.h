
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

#include "graph_io.h"
#include "mis/mis_config.h"
#include "mis/ils/ils.h"

struct instance {
  redu_vcc reduVCC;
  std::vector<reducer> reducer_stack;
  std::vector<instance> children;
};

class branch_and_reduce {
private:

  // std::vector<reducer> reducer_stack;

  std::string redu_type;
  std::string prune_type;
  std::string next_node_type;
  std::string enum_type;

  MISConfig config;

  NodeID min_deg_node(instance &inst);
  std::vector<std::vector<NodeID>> sorted_enumerate(instance &inst, NodeID x);

  void reduce(graph_access &G, instance &inst, reducer &R, unsigned int &num_fold_cliques, vertex_queue *queue);
  bool prune(instance &inst, unsigned int &curr_cover_size);
  std::vector<std::vector<NodeID>> enum_vertex(instance &inst, NodeID &v);
  NodeID nextNode(instance &inst);

  vertex_queue* construct_queue(graph_access &G, instance &inst,std::vector<NodeID> &clique);

  void buildCover(graph_access &G, instance *&inst);

  public:

    // redu_vcc reduVCC;
    instance root;

    unsigned int branch_count;
    std::vector<unsigned int> iso_degree;
    std::vector<unsigned int> dom_degree;
    unsigned int num_reductions;
    unsigned int num_attempts;

    // branch_and_reduce(graph_access &G);
    branch_and_reduce(graph_access &G, PartitionConfig &partition_config);
    virtual ~branch_and_reduce() {};

    void construct_run(PartitionConfig &partition_config);


    std::vector<std::vector<NodeID>> enumerate(instance &inst, NodeID v);
    void pivot_enumerator(instance &inst, std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);

    std::vector<std::vector<NodeID>> sorted_enumerate(instance &inst, NodeID x, std::vector<bool> &indset);

    void bandr( graph_access &G, instance &inst, unsigned int num_fold_cliques,
                vertex_queue *queue, PartitionConfig &partition_config, timer &t);


    void reduce_bnr( graph_access &G, instance &inst, unsigned int num_fold_cliques,
                 vertex_queue *queue, PartitionConfig &partition_config, timer &t);
    void branch_bnr( graph_access &G, instance &inst, unsigned int num_fold_cliques,
                 reducer &R, PartitionConfig &partition_config, timer &t);


    // void brute_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void reduMIS_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void small_degree_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void sort_enum_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);
    // void chalupa_status_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);
    // void cascading_red_bandr( graph_access &G, unsigned int num_fold_cliques, vertex_queue *queue,
                              // PartitionConfig &partition_config, timer &t);
    // void generate_mis_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t) {
      root.reduVCC.analyzeGraph(filename, G, t);
      std::cout << "Branches: " << branch_count << std::endl;
      std::cout << "ISO degree distrib. : [";
      for (unsigned int i : iso_degree) std::cout << i << ", ";
      std::cout << "]" << std::endl;
    };
};

#endif
