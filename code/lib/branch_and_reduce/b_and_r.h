
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

  std::vector<reducer> reducer_stack;

  std::string redu_type;
  std::string prune_type;
  std::string next_node_type;
  std::string enum_type;

  void reduce(graph_access &G, reducer &R, unsigned int &num_fold_cliques, vertex_queue *queue);
  bool prune(unsigned int &curr_cover_size);
  NodeID nextNode();
  NodeID min_deg_node();
  vertex_queue* construct_queue(graph_access &G, std::vector<std::vector<NodeID>> &clique);


  public:

    redu_vcc reduVCC;

    unsigned int branch_count;
    std::vector<unsigned int> iso_degree;
    std::vector<unsigned int> dom_degree;
    unsigned int num_reductions;
    unsigned int num_attempts;

    // branch_and_reduce(graph_access &G);
    branch_and_reduce(graph_access &G, PartitionConfig &partition_config);
    virtual ~branch_and_reduce() {};

    void construct_run(PartitionConfig &partition_config);


    std::vector<std::vector<NodeID>> enumerate(NodeID v);
    void pivot_enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);

    std::vector<std::vector<NodeID>> sorted_enumerate(NodeID x, std::vector<bool> &indset);

    void bandr( graph_access &G, unsigned int num_fold_cliques,
                vertex_queue *queue, PartitionConfig &partition_config, timer &t);



    // void brute_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void reduMIS_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void small_degree_bandr( graph_access &G, unsigned int num_fold_cliques);
    // void sort_enum_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);
    // void chalupa_status_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);
    // void cascading_red_bandr( graph_access &G, unsigned int num_fold_cliques, vertex_queue *queue,
                              // PartitionConfig &partition_config, timer &t);
    // void generate_mis_bandr( graph_access &G, unsigned int num_fold_cliques, PartitionConfig &partition_config, timer &t);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t) {reduVCC.analyzeGraph(filename, G, t);};
};

#endif
