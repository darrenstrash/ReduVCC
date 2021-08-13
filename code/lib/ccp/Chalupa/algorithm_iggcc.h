#ifndef ALGORITHM_IGGCC_H
#define ALGORITHM_IGGCC_H
#include "algorithm.h"
#include "timer.h"

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
//    long permutation_indset[MAX_VERTICES];
    long *permutation_indset;
//    long new_permutation_indset[MAX_VERTICES];
    long *new_permutation_indset;
//    long class_sizes[MAX_VERTICES];
    long *class_sizes;
//    long clique_sizes[MAX_VERTICES];
    long *clique_sizes;
//    long size_copies[MAX_VERTICES];
    long *size_copies;
//    bool changed_sizes[MAX_VERTICES];
    bool *changed_sizes;
//    long internal_positions[MAX_VERTICES];
    long *internal_positions;
//    long conflicting_vertices_list[MAX_VERTICES];
    long *conflicting_vertices_list;
//    long positioning_list[MAX_VERTICES];
    long *positioning_list;
//    long class_beginnings[MAX_VERTICES];
    long *class_beginnings;
//    long auxiliary_state[MAX_VERTICES];
    long *auxiliary_state;
//    long possibilities[MAX_VERTICES];
    long *possibilities;
//    long permutation[MAX_VERTICES];
    long *permutation;
//    long new_ordering[MAX_VERTICES];
    long *new_ordering;
//    vertex_list vertex_lists[MAX_VERTICES+1];
    vertex_list *vertex_lists;
//    bool have_adjacent_indset_member[MAX_VERTICES];
    bool *have_adjacent_indset_member;
//    long new_permutation[MAX_VERTICES];
    long *new_permutation;
//    long new_class_beginnings[MAX_VERTICES];
    long *new_class_beginnings;
    long generate_state_inverse_greedy(long permutation[]);
    long greedy_indset(graph G, long permutation_indset[]);
    void transform_representation();
//    bool free_indset_vertices[MAX_VERTICES];
    bool *free_indset_vertices;
//    long triangles_count[MAX_VERTICES];
    long *triangles_count;
    t_data *my_t_data;
public:
    algorithm_iggcc();
    ~algorithm_iggcc();
    bool iggcc_ccp(graph G, refer *result, refer *indset_size, refer *initial_indset, refer initial_indset_size, timer &total_timer, double &time_to_solution, int t_limit, int mis, std::size_t clique_cover_offset);
    bool gcc_ccp(graph G, refer *result, long *permutation);
    unsigned long long get_iterations();    
};

#endif // ALGORITHM_IGGCC_H
