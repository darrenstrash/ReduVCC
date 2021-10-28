
/******************************************************************************
 * reducer.h
 * * Preforms various reductions on reduVCC
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
    reducer();
    reducer(unsigned int iso_lim);
    virtual ~reducer() {};

    unsigned int num_reductions;
    unsigned int num_attempts;
    unsigned int num_cliques;
    unsigned int num_fold_cliques;


    // construct C from C'
    void unwindReductions(redu_vcc &reduVCC, double &time_to_solution);
    void unwindReductions(redu_vcc &reduVCC);
    // undo num reductions, constructing G from G'
    void undoReductions(redu_vcc &reduVCC);

    // brute appliciation of individual reductions
    void bruteISO(redu_vcc &reduVCC, std::vector<unsigned int> &iso_degree);
    void bruteD2(redu_vcc &reduVCC);
    void bruteTWIN(redu_vcc &reduVCC);
    void bruteDOM(redu_vcc &reduVCC, std::vector<unsigned int> &dom_degree);
    void bruteCROWN(redu_vcc &reduVCC);

    void exhaustive_reductions(redu_vcc &reduVCC,
                               std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &d2_degree);
    void cascading_reductions(redu_vcc &reduVCC, vertex_queue *queue,
                              std::vector<unsigned int> &iso_degree, std::vector<unsigned int> &d2_degree);

    std::size_t get_cover_size_offset() const;

};

#endif
