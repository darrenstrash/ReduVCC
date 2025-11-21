#include "algorithm_greedyindset.h"

algorithm_greedyindset::algorithm_greedyindset()
{
    have_adjacent_indset_member = new bool[MAX_VERTICES];
    Q = new refer[MAX_VERTICES];
    D = new double[MAX_VERTICES];
}

algorithm_greedyindset::~algorithm_greedyindset() {
    delete[] have_adjacent_indset_member;
    delete[] Q;
    delete[] D;
}

long algorithm_greedyindset::greedy_indset(graph G, refer result[])
{
    refer i,j;
    refer vertex;
    refer indset_size;

    this->G = G;
    createset();

    for (i=0;i<G->n;i++)
    {
        have_adjacent_indset_member[i] = false;
    }

    indset_size = 0;

    for (i=0;i<G->n;i++)
    {
        vertex = extractmin();
        // now we look, whether the next guy is not adjacent to some of the indset members
        if (! have_adjacent_indset_member[vertex])
        {
            result[indset_size] = vertex;
            indset_size++;
            for (j=0;j<G->V[vertex].edgecount;j++)
            {
                have_adjacent_indset_member[G->V[vertex].sibl[j]] = true;
            }
        }
    }

    return indset_size;
}

// vytvorenie prioritneho frontu
void algorithm_greedyindset::createset()
{
    refer i;
    n = G->n;
    for (i=0;i<n;i++)
    {
        Q[i] = i;
        // the priority = degree + random(0,1)
        D[i] = (double)(G->V[i].edgecount) + generator.random_double();
    }
    heapsize = n-1;
    buildheap();
}

// lavy potomok
refer algorithm_greedyindset::left(refer i)
{
    return 2*i+1;
}

// pravy potomok
refer algorithm_greedyindset::right(refer i)
{
    return 2*i+2;
}

// udrziavanie haldy
void algorithm_greedyindset::heapify(refer i)
{
    long l,r,largest;
    l = left(i);
    r = right(i);
    if (l <= heapsize && D[Q[l]] < D[Q[i]])
    {
        largest = l;
    }
    else
    {
        largest = i;
    }
    if (r <= heapsize && D[Q[r]] < D[Q[largest]])
    {
        largest = r;
    }
    if (largest != (long) i)
    {
        refer hlp;
        hlp = Q[i];
        Q[i] = Q[largest];
        Q[largest] = hlp;
        heapify(largest);
    }
}

// vytvorenie haldy
void algorithm_greedyindset::buildheap()
{
    long i;
    for (i=(long)(heapsize/2);i>=0;i--)
    {
        heapify(i);
    }
}

// extrakcia minima
refer algorithm_greedyindset::extractmin()
{
    refer min;
    if (heapsize < 0)
    {
        return 0;
    }
    min = Q[0];
    Q[0] = Q[heapsize];
    heapsize--;
    heapify(0);
    return min;
}

