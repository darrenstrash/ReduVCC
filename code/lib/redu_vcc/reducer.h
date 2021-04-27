
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

    unsigned int iso_limit;


  public:
    reducer() {};
    reducer(graph_access &G);
    reducer(graph_access &G, unsigned int iso_lim);
    virtual ~reducer() {};

    unsigned int num_reductions;
    unsigned int num_attempts;
    unsigned int num_cliques;
    unsigned int num_fold_cliques;


    // construct C from C'
    void unwindReductions(graph_access &G, redu_vcc *reduVCC);
    // undo num reductions, constructing G from G'
    void undoReductions(graph_access &G, redu_vcc &reduVCC);

    void bruteISO(graph_access &G, redu_vcc &reduVCC, std::vector<unsigned int> &iso_degree);
    void bruteD2(graph_access &G, redu_vcc &reduVCC);
    void bruteTWIN(graph_access &G, redu_vcc &reduVCC);
    void bruteDOM(graph_access &G, redu_vcc &reduVCC, std::vector<unsigned int> &dom_degree);
    void bruteCROWN(graph_access &G, redu_vcc &reduVCC);

    void exhaustive_reductions(graph_access &G, redu_vcc &reduVCC,
                               std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &d2_degree);
    void cascading_reductions(graph_access &G, redu_vcc &reduVCC, vertex_queue *queue,
                              std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &d2_degree);

};

#endif
