#include "common.h"

t_data my_t_data;

t_data *get_t_data()
{
    return &my_t_data;
}

void log(char message[])
{
    FILE *f;
    f = fopen("./log.txt","a");
    fprintf(f,"%s\n",message);
    fclose(f);
}

void swap(refer *a, refer *b)
{
        long pom;
        pom = *a;
        *a = *b;
        *b = pom;
}

long split(refer array[], long left, long right)
{
    long i = left;
    long j = right+1;
    long pivot = array[left];
    do
    {
        do
        {
            i++;
        }
        while (array[i] <= pivot && i <= right);
        do
        {
            j--;
        }
        while (array[j] > pivot);
        if (i < j)
        {
            swap(&array[i],&array[j]);
        }
    }
    while (i < j);
    swap(&array[left],&array[j]);
    return j;
}

void QuickSort(refer array[], long left, long right)
{
    long j;
    if (left < right)
    {
        j = split(array,left,right);
        QuickSort(array,left,j-1);
        QuickSort(array,j+1,right);
    }
}

