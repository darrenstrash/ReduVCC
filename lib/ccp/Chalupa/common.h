#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#define MAX_T_MAX 5100
#define MAX_INSTANCE 1
#define MAX_ALGORITHM 1

// datovy typ pre index vrchola
typedef unsigned int refer;
typedef struct {
    unsigned long t_cliques[MAX_T_MAX][6];
    unsigned long t_contributions[MAX_T_MAX];
} t_data;

t_data *get_t_data();
void log(char message[]);
void QuickSort(refer array[], long left, long right);

#endif
