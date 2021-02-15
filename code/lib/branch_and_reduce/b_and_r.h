
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef B_AND_R
#define B_AND_R

#include <algorithm>
#include "data_structure/graph_access.h"
#include "redu_vcc/reducer.h"

class branch_and_reduce {
private:

  reducer R;

  std::vector<bool> visited_nodes;

  public:
    branch_and_reduce(graph_access &G);
    virtual ~branch_and_reduce() {};

    std::vector<std::vector<NodeID>> enumerate(NodeID v);
    void enumerator(redu_vcc &reduVCC, std::vector<std::vector<NodeID>> &minimal_cliques,
                    std::vector<NodeID> &consider_nodes, std::vector<NodeID> &curr_clique, std::vector<NodeID> &excluded_nodes);

    void brute( graph_access &G);
    void branch( graph_access &G, unsigned int num_folded_cliques, NodeID curr_node);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
};

#endif
