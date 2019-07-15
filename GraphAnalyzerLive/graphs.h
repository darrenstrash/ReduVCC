#ifndef GRAPHS_H
#define GRAPHS_H
#include <stdio.h>
#include "common.h"
// maximalny pocet vrcholov grafu
////#define MAX_VERTICES 30050
#define MAX_VERTICES 5000000
////#define MAX_LABELS_CCP 40000
#define MAX_LABELS_CCP 6000000

/*
        ----------------
        DEKLARACIA GRAFU
        ----------------
*/

// vrchol grafu
typedef struct VERTEX
{
    refer edgecount;
    refer *sibl;
} vertex;

// samotny graf
typedef struct GRAPH_DATA
{
    refer n;
    unsigned long m;
    double density;
    vertex V[MAX_VERTICES];
} graph_data;
typedef graph_data *graph;

void input_graph(FILE *source);
// dealokujeme potrebnu pamat
void free_graph();
graph get_graph();
long get_problem();
bool are_adjacent(refer v, refer w);
void generate_graph_BA_model(unsigned long w, unsigned long n_max);

#endif // GRAPHS_H
