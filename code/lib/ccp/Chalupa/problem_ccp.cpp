#include "problem_ccp.h"

problem_ccp::problem_ccp()
{
}


unsigned long long problem_ccp::count_conflicts(graph G, refer *result)
{
    refer i,v;
    unsigned long long conflicts;

    class_sizes = new refer[G->n+1];

    conflicts = 0;

    for (v=0;v<G->n;v++)
    {
        class_sizes[v] = 0;
    }
    for (v=0;v<G->n;v++)
    {
        class_sizes[result[v]]++;
    }

    for (v=0;v<G->n;v++)
    {
        conflicts += class_sizes[result[v]]-1;

        for (i=0;i<G->V[v].edgecount;i++)
        {
            if (result[G->V[v].sibl[i]] == result[v])
            {
                conflicts--;
            }
        }
    }

    conflicts /= 2;

    return conflicts;
}

refer problem_ccp::find_conflicting_vertices(graph G, refer *result, refer *conflicting_vertices_list, refer *positioning_list)
{
    refer i,v;
    long f;
    unsigned long long conflicts;

    // computing the class sizes
    class_sizes = new refer[G->n+1];
    for (v=0;v<G->n;v++)
    {
        class_sizes[v] = 0;
    }
    for (v=0;v<G->n;v++)
    {
        class_sizes[result[v]]++;
    }

    conflicts = 0;
    for (v=0;v<G->n;v++)
    {
        f = class_sizes[result[v]]-1;

        for (i=0;i<G->V[v].edgecount;i++)
        {
            if (result[G->V[v].sibl[i]] == result[v])
            {
                f--;
            }
        }

        if (f > 0)
        {
            conflicting_vertices_list[conflicts] = v;
            positioning_list[v] = conflicts;
            conflicts++;
        }
    }

    delete(class_sizes);

    return conflicts;
}

refer problem_ccp::find_unlabeled_vertices(graph G, refer *result, refer *unlabeled_vertices_list)
{
    refer v,unlabeled;
    unlabeled = 0;
    for (v=0;v<G->n;v++)
    {
        if (result[v] <= 0)
        {
            unlabeled_vertices_list[unlabeled] = v;
            unlabeled++;
        }
    }
    return unlabeled;
}

void problem_ccp::label_unlabeled_vertices_randomly(graph G, refer *result, refer max_label)
{
    refer v;
    for (v=0;v<G->n;v++)
    {
        if (result[v] <= 0)
        {
            result[v] = generator.random(1,max_label);
        }
    }
}

refer problem_ccp::count_labels(graph G, refer *result)
{
    refer v,max_label;
    max_label = 0;
    for (v=0;v<G->n;v++)
    {
        //if (result[v] > G->n) {printf("WAA"); getchar();}
        if (max_label < result[v])
        {
            max_label = result[v];            
        }        
    }
    //getchar();
    return max_label;
}

