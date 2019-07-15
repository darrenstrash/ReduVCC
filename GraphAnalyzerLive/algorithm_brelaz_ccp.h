#ifndef ALGORITHM_BRELAZ_CCP_H
#define ALGORITHM_BRELAZ_CCP_H
#include "algorithm.h"
#include "graphs.h"
#include "vertex_label_function_refer.h"

class algorithm_brelaz_ccp : public algorithm
{
private:
    long long priorities[MAX_VERTICES];
    long long priorities_secondary[MAX_VERTICES];
    //short (*neighbor_color_matrix)[MAX_LABELS_CCP];
    vertex_label_function_refer *vertex_label_function;
    unsigned long class_sizes[MAX_VERTICES];
public:
    algorithm_brelaz_ccp(graph G, refer vertices);
    ~algorithm_brelaz_ccp();
    void brelaz_ccp(graph G, refer *result, refer max_labels);
};

#endif // ALGORITHM_BRELAZ_CCP_H
