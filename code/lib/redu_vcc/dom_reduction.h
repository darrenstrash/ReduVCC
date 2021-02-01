
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef DOM_RED
#define DOM_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "reduction.h"

class dom_reduction: public reduction {

  public:
    NodeID v;
    NodeID u;

    static bool validDOM(redu_vcc &reduVCC, NodeID &v, NodeID &u);
    static bool nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &a);

    void reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC);

};

#endif
