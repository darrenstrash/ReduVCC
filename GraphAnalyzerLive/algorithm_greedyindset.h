#ifndef ALGORITHM_GREEDYINDSET_H
#define ALGORITHM_GREEDYINDSET_H
#include "algorithm.h"
#include "graphs.h"

class algorithm_greedyindset : public algorithm
{
private:
    graph G;
    bool have_adjacent_indset_member[MAX_VERTICES];
    // halda
    refer Q[MAX_VERTICES];
    // D[i] is the priority value for i-th element of the heap
    double D[MAX_VERTICES];
    // velkost haldy
    refer heapsize;
    void createset();
    refer left(refer i);
    refer right(refer i);
    void heapify(refer i);
    void buildheap();
    refer extractmin();
public:
    algorithm_greedyindset();
    long greedy_indset(graph G, refer result[]);
};

#endif // ALGORITHM_GREEDYINDSET_H
