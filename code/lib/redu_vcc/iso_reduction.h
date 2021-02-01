
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef ISO_RED
#define ISO_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "reduction.h"

class iso_reduction: public reduction {

  public:
    NodeID v;

    static bool validISO(redu_vcc &reduVCC, NodeID &v);
    static bool validNeighbor(redu_vcc &reduVCC, NodeID &v, NodeID &u);

    void reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC) {};

};

#endif
