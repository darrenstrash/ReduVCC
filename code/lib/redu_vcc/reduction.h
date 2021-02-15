
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef REDU
#define REDU

#include <algorithm>
#include "data_structure/graph_access.h"
#include "redu_vcc.h"

class reduction {
  public:
    reduction() {};
    virtual ~reduction() {};

    // reduces the graph --> produce G'
    virtual void reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u ) = 0;
    // unfolds reductions --> produce C from C', maintain G'
    virtual void unfold(graph_access &G, redu_vcc &reduVCC) = 0;
    // undoes reduction --> produces G from G'
    virtual void unreduce(graph_access &G, redu_vcc &reduVCC) = 0;

    bool isSubset(redu_vcc &reduVCC, std::vector<NodeID> &A, std::vector<NodeID> &B);
    void merge_neighborhoods(redu_vcc &reduVCC, NodeID &a, NodeID &b);

    static bool uncrossedSets(redu_vcc &reduVCC, NodeID &a, NodeID &b);

};

#endif
