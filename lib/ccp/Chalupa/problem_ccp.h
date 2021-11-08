#ifndef PROBLEM_CCP_H
#define PROBLEM_CCP_H
#include "problem.h"
#include "graphs.h"

class problem_ccp : public problem
{
private:
    refer *class_sizes;
public:
    problem_ccp();
    unsigned long long count_conflicts(graph G, refer *result);
    refer find_conflicting_vertices(graph G, refer *result, refer *conflicting_vertices_list, refer *positioning_list);
    refer find_unlabeled_vertices(graph G, refer *result, refer *unlabeled_vertices_list);
    void label_unlabeled_vertices_randomly(graph G, refer *result, refer max_label);
    refer count_labels(graph G, refer *result);
};

#endif // PROBLEM_CCP_H
