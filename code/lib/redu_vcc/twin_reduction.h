
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/

#ifndef TWIN_RED
#define TWIN_RED

#include "data_structure/graph_access.h"
#include "redu_vcc.h"
#include "reduction.h"

class twin_reduction: public reduction {

  private:
    std::vector<NodeID> edge_nodes;
    NodeID nonedge_node;

    std::vector<NodeID> N_w; // neighborhood of w before fold, excluding v
    std::vector<NodeID> N_x; // neighborhood of x before fold, excluding v

    std::vector<NodeID> disjoint; // vertices in N_w and N_x not in N_y

    bool remove_type;
    std::vector<NodeID> clique1;
    std::vector<NodeID> clique2;

    bool removeType (redu_vcc &reduce);
    void removeTWIN (redu_vcc &reduVCC);
    void foldTWIN(redu_vcc &reduVCC);

    void unfoldTWIN(redu_vcc &reduVCC, std::vector<NodeID> &partial_clique,
                    unsigned int &clique_id, NodeID &a, NodeID &b, NodeID &c);

  public:
    NodeID v;
    NodeID u;
    NodeID w;
    NodeID x;
    NodeID y;


    static bool validTWIN(redu_vcc &reduVCC, NodeID &v, NodeID &u);
    static void assignNodes(redu_vcc &reduVCC, NodeID &v, NodeID &w, NodeID &x, NodeID &y);
    static bool twinFound( redu_vcc &reduVCC,  NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y);
    static bool validNeighbors(redu_vcc &reduVCC, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y);

    void reduce(graph_access &G, redu_vcc &reduVCC, NodeID &node_v, NodeID &node_u );
    void unreduce(graph_access &G, redu_vcc &reduVCC);
    void unfold(graph_access &G, redu_vcc &reduVCC);

};

#endif
