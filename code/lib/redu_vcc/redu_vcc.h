
/******************************************************************************
 * redu_vcc.h
 * *
 *
 *****************************************************************************/

#ifndef REDU_VCC
#define REDU_VCC

#include "data_structure/graph_access.h"


#define for_adjList(reduVCC, v, n) { for (NodeID n : reduVCC.adj_list[v]) { if (!reduVCC.node_status[n]) continue;

class redu_vcc {
  private:
    // subgraph mapping
    std::vector<int> old_to_new_map;
    std::vector<NodeID> new_to_old_map;
    void assignMaps(graph_access &G);

    // generate adj_list from G
    void generateAdjList(graph_access &G);
    // checks if clique is in original graph G
    bool cliqueInG(graph_access &G, std::vector<NodeID> &clique);




  public:
    std::vector<std::vector<NodeID>> adj_list;
    std::vector<bool> node_status;  // marks nodes status in G
    std::vector<bool> fold_node;   // marks nodes removed in a fold
    unsigned int remaining_nodes;

    std::vector<std::vector<NodeID>> clique_cover;

    unsigned int next_cliqueID; // next cliqueID while reducing
    unsigned int next_solvecliqueID; // next cliqueID while solving
    std::vector<unsigned int> node_clique; // node to cliqueID mapping
    std::vector<unsigned int> solve_node_clique; // node to cliqueID mapping for solving

    std::vector<std::vector<int>> kernel_adj_list; // subgraph adj_list
    unsigned long kernel_edges; // num subgraph edges

    std::vector<bool> scratch1;
    std::vector<bool> scratch2;

    redu_vcc() {};
    redu_vcc(graph_access &G);
    virtual ~redu_vcc() {};

    void build_cover(graph_access &G);
    void validateCover(graph_access &G);

    void buildKernel(graph_access &G);
    void addKernelCliques(std::vector<std::vector<int>> &clique_set);
    void addCrownCliques(std::vector<std::vector<NodeID>> &crown_cliques, std::vector<std::vector<int>> &clique_set);


    unsigned int adj_size(NodeID v);
    std::vector<NodeID> curr_adj_list(NodeID v);

    void removeVertex(NodeID v);
    void addVertex(NodeID v);
    void removeVertexSet(std::vector<NodeID> &S);
    void addVertexSet(std::vector<NodeID> &S);

    void addClique(std::vector<NodeID> &clique);
    void addCliqueToCover(std::vector<NodeID> &clique);
    // std::vector<NodeID> pop_clique();
    void pop_clique(std::vector<NodeID> &clique);


    unsigned int getCliqueID(NodeID &v);
    std::vector<NodeID> getClique(NodeID &v);
    void replaceClique(unsigned int cliqueID, std::vector<NodeID> new_clique);

    // void clearScratch(std::vector<bool> &scratch);

    void printAdjList();
    void printAdjList(NodeID v);
    void printNeighborhood(NodeID v);
    void printVectorSet(std::vector<NodeID> S);

};

#endif
