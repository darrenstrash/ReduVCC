
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef UNCONF_RED
#define UNCONF_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "vertex_queue.h"
#include "reduction.h"

class unconfined_reduction: public reduction {

  public:
    NodeID v;
    NodeID u;

    std::vector<NodeID> clique;

    static bool validUNCONFINED(redu_vcc &reduVCC, NodeID &v, NodeID &u);
    // static bool nodeDominates(redu_vcc &reduVCC, NodeID &v, NodeID &a);

    void reduce(redu_vcc &reduVCC,
                NodeID &node_v, NodeID &node_u );
    void reduce( redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unreduce(redu_vcc &reduVCC);
    void unfold( redu_vcc &reduVCC);

};

#endif
