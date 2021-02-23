
/******************************************************************************
 * reducer.h
 * *
 *
 *****************************************************************************/
 #ifndef REDUCER
 #define REDUCER


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

    reducer() {};
    virtual ~reducer() {};

    void init(graph_access &G);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);

    // map node_clique to clique_cover
    void buildCover(graph_access &G);
    // solve remaining_nodes with Chalupa Solver
    void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);
    // construct C from C'
    void unwindReductions(graph_access &G);
    // undo num reductions, constructing G from G'
    void undoReductions(graph_access &G, unsigned int num);

    std::vector<unsigned int> bruteISO(graph_access &G);
    std::vector<unsigned int> bruteD2(graph_access &G);
    std::vector<unsigned int> bruteTWIN(graph_access &G);
    std::vector<unsigned int> bruteDOM(graph_access &G);
    std::vector<unsigned int> bruteCROWN(graph_access &G);

};

#endif
