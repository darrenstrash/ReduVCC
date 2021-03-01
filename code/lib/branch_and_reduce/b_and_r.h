
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef B_AND_R
#define B_AND_R

#include <algorithm>
#include <functional>

#include "timer.h"
#include "data_structure/graph_access.h"
#include "redu_vcc/reducer.h"

class branch_and_reduce {
private:

  reducer R;

  std::vector<bool> visited_nodes;

  public:
    unsigned int mis;
    std::vector<bool> node_mis;

    branch_and_reduce(graph_access &G);
    virtual ~branch_and_reduce() {};

    void solveKernel(graph_access &G, PartitionConfig partition_config, timer &t);
    void getMIS(std::string file);

    unsigned int remainingMIS(redu_vcc &reduVCC);
    unsigned int overlapMIS(std::vector<NodeID> &clique);

    unsigned int exhaustive_reductions(graph_access &G, unsigned int &num_folded_cliques, unsigned int &curr_mis);

    std::vector<std::vector<NodeID>> enumerate(NodeID v);
    void enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);
    void pivot_enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);

    void brute( graph_access &G);
    void prune(graph_access &G);
    void min_degree_prune(graph_access &G);
    void upper_bound_prune(graph_access &G, PartitionConfig &partition_config, timer &s);

    void branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node);
    void prune_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node);
    void small_deg_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node);
    void bounding_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis, NodeID curr_node);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
};

#endif
