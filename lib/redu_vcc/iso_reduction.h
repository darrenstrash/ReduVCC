
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

    void reduce( redu_vcc &reduVCC,
                NodeID &node_v, NodeID &node_u );
    void reduce( redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unfold( redu_vcc &reduVCC);
    void unreduce( redu_vcc &reduVCC);

};

#endif
