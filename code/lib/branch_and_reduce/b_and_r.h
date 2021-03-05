
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

  public:

    redu_vcc reduVCC;

    unsigned int mis;
    std::vector<bool> node_mis;

    std::vector<reducer> reducer_stack;
    unsigned int num_reductions;

    branch_and_reduce(graph_access &G);
    virtual ~branch_and_reduce() {};

    void getMIS(std::string file);

    unsigned int remainingMIS();
    unsigned int overlapMIS(std::vector<NodeID> &clique);

    // unsigned int exhaustive_reductions(graph_access &G, unsigned int &num_folded_cliques, unsigned int &curr_mis);

    std::vector<std::vector<NodeID>> enumerate(NodeID v);
    void enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);
    void pivot_enumerator(std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);

    // void brute( graph_access &G);
    // void prune(graph_access &G);
    // void min_degree_prune(graph_access &G);
    //
    void branch( graph_access &G, unsigned int num_folded_cliques);
    void prune_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis);
    void small_deg_branch( graph_access &G, unsigned int num_folded_cliques, unsigned int curr_mis);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t) {reduVCC.analyzeGraph(filename, G, t);};
};

#endif
