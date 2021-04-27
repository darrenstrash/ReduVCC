
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef DOM_RED
#define DOM_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "vertex_queue.h"
#include "reduction.h"

class dom_reduction: public reduction {

  public:
    NodeID v;
    NodeID u;

    static bool validDOM(redu_vcc &reduVCC, NodeID &v, NodeID &u);
    static bool nodeDominates(redu_vcc &reduVCC, NodeID &v, NodeID &a);

    void reduce(graph_access &G, redu_vcc &reduVCC,
                NodeID &node_v, NodeID &node_u );
    void reduce(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC);
    void unfold(graph_access &G, redu_vcc* reduVCC);

};

#endif
