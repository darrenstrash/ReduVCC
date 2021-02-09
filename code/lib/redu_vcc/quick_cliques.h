
/******************************************************************************
 * quick_cliques.h
 * *
 *
 *****************************************************************************/

#ifndef QUICK_CLIQUES
#define QUICK_CLIQUES

#include "data_structure/graph_access.h"

class quick_cliques {
  private:
    std::vector<int> old_to_new_map;
    std::vector<NodeID> new_to_old_map;

    void assignMaps(graph_access &G);



  public:
    std::vector<std::list<int>> int_adj_list;

    quick_cliques() {};
    virtual ~quick_cliques() {};

    void buildIntAdjList(graph_access &G);

};

#endif
