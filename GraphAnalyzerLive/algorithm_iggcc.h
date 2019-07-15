#ifndef ALGORITHM_IGGCC_H
#define ALGORITHM_IGGCC_H
#include "algorithm.h"

class vertex_list
{
public:
    bool occ;
    long vertex;
    vertex_list *next;
    vertex_list *prev;
    vertex_list *end;
};

class algorithm_iggcc : public algorithm
{
private:
    graph G;
    refer *result;
    refer colors;
    long conflicts;
    long max_size;
    long current_class;
    long fitness_indset;
    long new_fitness_indset;
    long permutation_indset[MAX_VERTICES];
    long new_permutation_indset[MAX_VERTICES];
    long class_sizes[MAX_VERTICES];
    long clique_sizes[MAX_VERTICES];
    long size_copies[MAX_VERTICES];
    bool changed_sizes[MAX_VERTICES];
    long internal_positions[MAX_VERTICES];
    long conflicting_vertices_list[MAX_VERTICES];
    long positioning_list[MAX_VERTICES];
    long class_beginnings[MAX_VERTICES];
    long auxiliary_state[MAX_VERTICES];
    long possibilities[MAX_VERTICES];
    long permutation[MAX_VERTICES];
    long new_ordering[MAX_VERTICES];
    vertex_list vertex_lists[MAX_VERTICES+1];
    bool have_adjacent_indset_member[MAX_VERTICES];
    long new_permutation[MAX_VERTICES];
    long new_class_beginnings[MAX_VERTICES];
    long generate_state_inverse_greedy(long permutation[]);
    long greedy_indset(graph G, long permutation_indset[]);
    void transform_representation();
    bool free_indset_vertices[MAX_VERTICES];
    long triangles_count[MAX_VERTICES];
    t_data *my_t_data;
public:
    algorithm_iggcc();
    bool iggcc_ccp(graph G, refer *result, refer *indset_size, refer *initial_indset, refer initial_indset_size);
    bool gcc_ccp(graph G, refer *result, long *permutation);
    unsigned long long get_iterations();    
};

#endif // ALGORITHM_IGGCC_H
