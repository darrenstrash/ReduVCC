
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
    reducer() {};
    reducer(graph_access &G);
    virtual ~reducer() {};

    unsigned int num_reductions;
    unsigned int num_cliques;
    unsigned int num_fold_cliques;


    // construct C from C'
    void unwindReductions(graph_access &G, redu_vcc &reduVCC);
    // undo num reductions, constructing G from G'
    void undoReductions(graph_access &G, redu_vcc &reduVCC);

    void bruteISO(graph_access &G, redu_vcc &reduVCC);
    void bruteD2(graph_access &G, redu_vcc &reduVCC);
    void bruteTWIN(graph_access &G, redu_vcc &reduVCC);
    void bruteDOM(graph_access &G, redu_vcc &reduVCC);
    void bruteCROWN(graph_access &G, redu_vcc &reduVCC);

    void exhaustive_reductions(graph_access &G, redu_vcc &reduVCC);
    void cascading_reductions(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue);

};

#endif
