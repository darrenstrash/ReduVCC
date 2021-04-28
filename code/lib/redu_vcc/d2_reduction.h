
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef D2_RED
#define D2_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "vertex_queue.h"
#include "reduction.h"

class d2_reduction: public reduction {

  private:
    void foldD2(redu_vcc &reduVCC);

  public:
    NodeID v;
    NodeID u;
    NodeID w;

    std::vector<NodeID> N_u; // neighborhood of u before fold, excluding v
    std::vector<NodeID> disjoint; // N_u \cap N_w

    static bool validD2(redu_vcc &reduVCC, NodeID &v);
    static void assignNodes(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID &w);
    static bool isTriangle(redu_vcc &reduVCC, NodeID &u, NodeID &w);

    void reduce(graph_access &G, redu_vcc &reduVCC,
                NodeID &node_v, NodeID &node_u );
    void reduce(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC);
    void unfold(graph_access &G, redu_vcc &reduVCC);

};

#endif
