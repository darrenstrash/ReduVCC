
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/
 #ifndef REDUCER
 #define REDUCER

// #include "quick-cliques/Tools.h"
// #include "quick-cliques/TomitaAlgorithm.h"
// #include "quick-cliques/AdjacencyListAlgorithm.h"
// #include "quick-cliques/HybridAlgorithm.h"
// #include "quick-cliques/DegeneracyAlgorithm.h"
// #include "quick-cliques/CliqueTools.h"


#include "partition/partition_config.h"
#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_vcc.h"
#include "reduction.h"
#include "iso_reduction.h"
#include "d2_reduction.h"
#include "twin_reduction.h"
#include "dom_reduction.h"
#include "crown_reduction.h"

class reducer {
  private:
    std::vector<reduction*> reduction_stack;


  public:
    redu_vcc reduVCC;

    reducer(graph_access &G);
    virtual ~reducer() {};

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);

    void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);
    void unwindReductions(graph_access &G);

    void bruteISO(graph_access &G);
    void bruteD2(graph_access &G);
    void bruteTWIN(graph_access &G);
    void bruteDOM(graph_access &G);
    void bruteCROWN(graph_access &G);

    std::vector<std::vector<NodeID>> enumerate(NodeID v);
    void enumerator(std::vector<NodeID> &N_v,
                    unsigned int curr_i,
                    std::vector<NodeID> curr_clique,
                    std::vector<std::vector<NodeID>> &enum_cliques);

  void branch_and_bound();
  void branch( std::vector<std::vector<NodeID>> &clique_cover,
                        std::vector<std::vector<NodeID>> &curr_cover,
                        NodeID curr_node);

};

#endif
