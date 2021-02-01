
/******************************************************************************
 * redu_vcc.h
 * *
 *
 *****************************************************************************/

#ifndef REDU_VCC
#define REDU_VCC

#include "data_structure/graph_access.h"

class redu_vcc {
  private:
    std::vector<int> old_to_new_map;
    std::vector<NodeID> new_to_old_map;

    void generateAdjList(graph_access &G);
    bool cliqueInG(graph_access &G, std::vector<NodeID> &clique);

    void assignMaps(graph_access &G);



  public:
    std::vector<std::vector<NodeID>> adj_list;
    std::vector<bool> node_status;
    unsigned int remaining_nodes;

    std::vector<std::vector<NodeID>> clique_cover;
    std::vector<unsigned int> node_clique;

    std::vector<std::vector<int>> kernel_adj_list;
    unsigned long kernel_edges;

    std::vector<bool> scratch1;
    std::vector<bool> scratch2;

    redu_vcc() {};
    virtual ~redu_vcc() {};

    void build(graph_access &G);
    void validateCover(graph_access &G);

    void buildKernel(graph_access &G);
    void addKernelCliques(std::vector<std::vector<int>> &clique_set);

    void removeVertex(NodeID v);

    void addClique(std::vector<NodeID> &clique);
    void removeVertexSet(std::vector<NodeID> &S);

    unsigned int getCliqueID(NodeID &v);
    std::vector<NodeID> getClique(NodeID &v);
    void replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique);

    void clearScratch(std::vector<bool> &scratch);

    void printAdjList();
    void printAdjList(NodeID v);
    void printNeighborhood(NodeID v);
    void printVectorSet(std::vector<NodeID> S);

};

#endif
