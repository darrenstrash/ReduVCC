#ifndef STATISTICS_H
#define STATISTICS_H
#include <math.h>
#include "graphs.h"

class statistics
{
public:
    statistics();
    static refer min_degree(graph G);
    static refer max_degree(graph G);
    static double average_degree(graph G);
    static double degree_stdev(graph G);
    static refer degree_distribution(graph G, long *target);
    static unsigned long triangles(graph G);
};

#endif // STATISTICS_H
