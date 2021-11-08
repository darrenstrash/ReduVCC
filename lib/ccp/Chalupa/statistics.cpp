#include "statistics.h"

statistics::statistics()
{
}

refer statistics::min_degree(graph G)
{
    long i;
    refer min_deg;

    min_deg = G->n;
    for (i=0;i<G->n;i++)
    {
        if (min_deg > G->V[i].edgecount)
        {
            min_deg = G->V[i].edgecount;
        }
    }

    return min_deg;
}

refer statistics::max_degree(graph G)
{
    long i;
    refer max_deg;

    max_deg = 0;
    for (i=0;i<G->n;i++)
    {
        if (max_deg < G->V[i].edgecount)
        {
            max_deg = G->V[i].edgecount;
        }
    }

    return max_deg;
}

double statistics::average_degree(graph G)
{
    long i;
    double avg_deg;

    avg_deg = 0.0;
    for (i=0;i<G->n;i++)
    {
        avg_deg += (double)(G->V[i].edgecount);
    }
    avg_deg /= G->n;

    return avg_deg;
}

double statistics::degree_stdev(graph G)
{
    long i;
    double avg_deg,stdev_deg;

    avg_deg = 0.0;
    for (i=0;i<G->n;i++)
    {
        avg_deg += (double)(G->V[i].edgecount);
    }
    avg_deg /= G->n;

    stdev_deg = 0.0;
    for (i=0;i<G->n;i++)
    {
        stdev_deg += ((double)(G->V[i].edgecount)-avg_deg)*((double)(G->V[i].edgecount)-avg_deg);
    }
    stdev_deg /= G->n;
    stdev_deg = sqrt(stdev_deg);

    return stdev_deg;
}

refer statistics::degree_distribution(graph G, long *target)
{
    // degree distr.
    long i, maxdeg = 0;

    for (i=0;i<G->n;i++)
    {
        target[i] = 0;
    }
    for (i=0;i<G->n;i++)
    {
        if (maxdeg < G->V[i].edgecount)
        {
            maxdeg = G->V[i].edgecount;
        }
        target[G->V[i].edgecount] += 1;
    }

    return maxdeg;
}

unsigned long statistics::triangles(graph G)
{
    long i,j,k,l;
    unsigned long count = 0;

    for (i=0;i<G->n;i++)
    {
        for (j=0;j<G->V[i].edgecount;j++)
        {
            for (k=0;k<G->V[G->V[i].sibl[j]].edgecount;k++)
            {
                for (l=0;l<G->V[G->V[G->V[i].sibl[j]].sibl[k]].edgecount;l++)
                {
                    if (i == G->V[G->V[G->V[i].sibl[j]].sibl[k]].sibl[l])
                    {
                        count++;
                    }
                }
            }
        }
    }

    return count / 2;
}
