#ifndef ALGORITHM_H
#define ALGORITHM_H
#include "random_generator.h"
#include "graphs.h"

class algorithm
{
protected:
    random_generator generator;
public:
    algorithm();
    virtual bool gcc_ccp(graph G, refer *result, long *permutation) { return false; }
};

#endif // ALGORITHM_H
