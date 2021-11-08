
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef CROWN_RED
#define CROWN_RED

#include "mis/kernel/branch_and_reduce_algorithm.h"
#include "data_structure/graph_access.h"
#include "vertex_queue.h"
#include "redu_vcc.h"
#include "reduction.h"

class crown_reduction: public reduction {

  public:
    std::vector<std::vector<NodeID>> crown_cliques;

    void reduce(redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u );
    void reduce(redu_vcc &reduVCC, vertex_queue *queue,
                NodeID &node_v, NodeID &node_u );
    void unreduce(redu_vcc &reduVCC);
    void unfold(redu_vcc &reduVCC) {};

};

#endif
