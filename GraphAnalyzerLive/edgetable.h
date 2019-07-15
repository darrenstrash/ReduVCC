#ifndef EDGETABLE_H
#define EDGETABLE_H
#include <stdio.h>
#include <stdlib.h>
#include "graphs.h"
#define MAX_HASH 500000

class hashset
{
public:
    bool occ;
    refer v;
    refer w;
    hashset *next;
};

class edgetable
{
private:
    hashset H[MAX_HASH];
public:
    edgetable();
    // hashovacia funkcia
    unsigned long long hash(refer v, refer w);
    bool isin(refer v, refer w);
    void insert(refer v, refer w);
    void clear();
};

#endif // EDGETABLE_H
