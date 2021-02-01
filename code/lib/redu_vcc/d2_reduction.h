
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef D2_RED
#define D2_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "reduction.h"

class d2_reduction: public reduction {

  private:
    void foldD2(redu_vcc &reduVCC);

  public:
    NodeID v;
    NodeID u;
    NodeID w;

    std::vector<NodeID> N_u; // neighborhood of u before fold, excluding v

    static bool validD2(redu_vcc &reduVCC, NodeID &v);
    static void assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w);
    static bool isTriangle(std::vector<std::vector<NodeID>> &adj_list, NodeID &u, NodeID &w);
    static bool validNeighbors(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID & w);

    void reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC);

};

#endif
