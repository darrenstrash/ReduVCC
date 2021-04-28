
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef ISO_RED
#define ISO_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "vertex_queue.h"
#include "reduction.h"

class iso_reduction: public reduction {

  public:
    NodeID v;
    std::vector <NodeID> clique;

    static bool validISO(redu_vcc &reduVCC, unsigned int& deg_limit, NodeID &v);
    static bool validNeighbor(redu_vcc &reduVCC, NodeID &v, NodeID &u);

    void reduce(graph_access &G, redu_vcc &reduVCC,
                NodeID &node_v, NodeID &node_u );
    void reduce(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unfold(graph_access &G, redu_vcc &reduVCC);
    void unreduce(graph_access &G, redu_vcc &reduVCC);

};

#endif
