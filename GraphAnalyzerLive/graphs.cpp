/*
        ===================
        IMPLEMENTACIA GRAFU
        ===================

        G = [V,E]
        |V| = n
        |E| = m
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "graphs.h"
#include "edgetable.h"
#include "random_generator.h"

#include <iostream>

graph G;
edgetable *edges;

/*
  ---------------
  NACITANIE GRAFU
  ---------------
*/

// pomocne polia pre jednotlive hrany grafu (aby stacil jeden malloc)
refer *E1;
refer *E2;
// pomocne pole obsahujuce stupne vrcholov
refer VEc[MAX_VERTICES];
unsigned long roulette[MAX_VERTICES];
unsigned long roulette_sum;
random_generator generator;

graph get_graph()
{
    return G;
}

// vstup grafu
void input_graph(FILE *source)
{
    char ch;
    unsigned long i, j;
    refer v,w,n;
    long long cnt = 0;

    // precitame komentare
    while ((ch = fgetc(source)) == 'c')
    {
        while ((ch = fgetc(source)) != '\n');
    }

    // overridneme retazec "p edge" alebo "p team"
    for (i=0;i<1;i++)
    {
        fgetc(source);
    }
    for (i=2;i<6;i++)
    {
        fgetc(source);
    }
    
    std::cout << "reads p edge" << std::endl;

    // alokujeme priestor pre graf
    G = (graph) malloc (sizeof(graph_data));

    // nacitame pocty vrcholov a hran
    fscanf(source,"%u",&G->n); fscanf(source,"%lu",&G->m);
    n = G->n;

    // allocation of the auxiliary space
    E1 = new refer[2*G->m+1];
    E2 = new refer[2*G->m+1];
    edges = new edgetable();

    // nacitame novy riadok
    fgetc(source);

    edges->clear();

    std::cout << "here" << std::endl;
    // nacitame hrany a vlozime ich do pomocneho pola
    long real_m = 0;
    for (i=0;i<G->m;i++)
    {
        // nacitame novy riadok a pismeno e
        fgetc(source);
        fgetc(source);
        if (fscanf(source,"%u",&v) == EOF)
        {
            break;
        }
        if (fscanf(source,"%u",&w) == EOF)
        {
            break;
        }
        std::cout << "read edge " << v << ", " << w << std::endl;
        // akceptujeme iba jeden "smer"
        v--; w--;
        //if (v >= 50000 || w >= 50000) {printf("fuj!(%ld)",v); getchar();}
        if ((! are_adjacent(v,w) || ! are_adjacent(w,v)) && v != w)
        {
            if (! are_adjacent(v,w) && ! are_adjacent(w,v))
            {
                real_m++;
            }
            if (! are_adjacent(v,w))
            {
                E1[cnt] = v; E2[cnt] = w;
                VEc[v]++; cnt++;
                edges->insert(v,w);
            }
            if (! are_adjacent(w,v))
            {
                E1[cnt] = w; E2[cnt] = v;
                VEc[w]++; cnt++;
                edges->insert(w,v);
            }
        }
        fgetc(source);
    }
    
    std::cout << "here2" << std::endl;
    //G->n = n;

    // alokujeme potrebnu pamat
    for (i=0;i<G->n;i++)
    {
        G->V[i].edgecount = 0;
        //printf("%u/%ld,",G->n,VEc[i]);
        G->V[i].sibl = (refer *) malloc (VEc[i]*sizeof(refer));
    }
    // vlozime potrebne udaje
    std::cout << "here3" << std::endl;
    for (j=0;j<cnt;j++)
    {
        v = E1[j]; w = E2[j];
        G->V[v].edgecount++;
        G->V[v].sibl[G->V[v].edgecount-1] = w;
    }
    G->m = real_m;
    G->density = ((double)(G->m))/((double)(G->n*(G->n-1)/2));

    std::cout << "here4" << std::endl;
    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
    }
    std::cout << "here5" << std::endl;
    delete(E1);
    delete(E2);
    edges->clear();
    delete(edges);
}

long blee = 0;
long blee_count = 0;

// generovanie grafu podla BA modelu
void generate_graph_BA_model(unsigned long w, unsigned long n_max)
{
    unsigned long i;
    unsigned long j,q,r;
    unsigned long triangles;
    refer triangles_involved[MAX_VERTICES];

    srand((unsigned long long) time(NULL));

    // alokujeme priestor pre graf
    G = (graph) malloc (sizeof(graph_data));

    edges = new edgetable();
    edges->clear();

    // inicializujeme zakladny graf
    G->n = w;
    G->m = w-1;
    for (i=0;i<n_max;i++)
    {
        G->V[i].sibl = (refer *) malloc (1000*sizeof(refer));
    }
    for (i=0;i<w;i++)
    {
        if (i == 0 || i == w - 1)
        {
            G->V[i].edgecount = 1;
            if (i == 0)
            {
                G->V[i].sibl[0] = i+1;
            }
            if (i == w - 1)
            {
                G->V[i].sibl[0] = i-1;
            }
        }
        else
        {
            G->V[i].edgecount = 2;
            G->V[i].sibl[0] = i-1;
            G->V[i].sibl[1] = i+1;
        }

    }
    // zly, zakerny zakladny graf
    /*G->n = 10;
    G->m = 17;
    for (i=0;i<n_max;i++)
    {
        G->V[i].sibl = (refer *) malloc (1000*sizeof(refer));
    }
    q = 0; G->V[0].edgecount = 4; G->V[q].sibl[0] = 1; G->V[q].sibl[1] = 2; G->V[q].sibl[2] = 3; G->V[q].sibl[3] = 5;
    q = 1; G->V[1].edgecount = 3; G->V[q].sibl[0] = 0; G->V[q].sibl[1] = 6; G->V[q].sibl[2] = 8;
    q = 2; G->V[2].edgecount = 4; G->V[q].sibl[0] = 0; G->V[q].sibl[1] = 3; G->V[q].sibl[2] = 6; G->V[q].sibl[3] = 7;
    q = 3; G->V[3].edgecount = 4; G->V[q].sibl[0] = 0; G->V[q].sibl[1] = 2; G->V[q].sibl[2] = 4; G->V[q].sibl[3] = 6;
    q = 4; G->V[4].edgecount = 3; G->V[q].sibl[0] = 3; G->V[q].sibl[1] = 5; G->V[q].sibl[2] = 7;
    q = 5; G->V[5].edgecount = 4; G->V[q].sibl[0] = 0; G->V[q].sibl[1] = 4; G->V[q].sibl[2] = 6; G->V[q].sibl[3] = 9;
    q = 6; G->V[6].edgecount = 5; G->V[q].sibl[0] = 1; G->V[q].sibl[1] = 2; G->V[q].sibl[2] = 3; G->V[q].sibl[3] = 5; G->V[q].sibl[4] = 9;
    q = 7; G->V[7].edgecount = 2; G->V[q].sibl[0] = 2; G->V[q].sibl[1] = 4;
    q = 8; G->V[8].edgecount = 3; G->V[q].sibl[0] = 1; G->V[q].sibl[1] = 9;
    q = 9; G->V[9].edgecount = 2; G->V[q].sibl[0] = 5; G->V[q].sibl[1] = 6; G->V[q].sibl[2] = 8;
    G->density = ((double)(G->m))/((double)(G->n*(G->n-1)/2));*/

    triangles = 0;

    refer n = G->n;
    // a teraz ideme na samotny BA model
    for (i=0;i<n_max;i++)
    {
        triangles_involved[i] = 0;
    }
    for (i=n;i<n_max;i++)
    {
        G->V[i].edgecount = w;
        // roulette wheel preparation
        roulette[0] = 0;
        roulette_sum = 0;
        for (j=1;j<i;j++)
        {
            roulette[j] = roulette[j-1] + G->V[j-1].edgecount;
            roulette_sum += G->V[j-1].edgecount;
        }
        roulette_sum += G->V[i-1].edgecount;
        // roulette wheel selection of new vertices
        for (j=0;j<w;j++)
        {
            do
            {
                r = generator.random(0,roulette_sum-1);
                q = 0;
                while (! (roulette[q] <= r && (q >= i-1 || roulette[q+1] > r)))
                {
                    q++;
                }
            }
            while (are_adjacent(i,q) || are_adjacent(q,i));
            G->V[i].sibl[j] = q;
            G->V[q].sibl[G->V[q].edgecount] = i;
            G->V[q].edgecount++;
            edges->insert(i,q);
            edges->insert(q,i);
        }
        G->n++;
        G->m += w;
        if (are_adjacent(G->V[G->n-1].sibl[0],G->V[G->n-1].sibl[1]))
        {
            triangles++;
            triangles_involved[G->n-1]++;
            triangles_involved[G->V[G->n-1].sibl[0]]++;
            triangles_involved[G->V[G->n-1].sibl[1]]++;
        }
        //printf("%ld,",triangles);
    }
    //getchar();
    /*for (i=0;i<G->n;i++)
    {
        if (triangles_involved[i] > 0)
        {
            printf("%u",triangles_involved[i]);
            printf("(%u),",G->V[i].edgecount);
        }
    }*/
    printf("%ld,",triangles);
    blee += triangles;
    blee_count++;
    printf("(%ld)%0.2lf,",triangles,(double)(blee)/(double)(blee_count));
    //getchar();

    for (i=0;i<G->n;i++)
    {
        QuickSort(G->V[i].sibl,0,G->V[i].edgecount-1);
        //printf("%ld->",i);
        /*for (j=0;j<G->V[i].edgecount;j++)
        {
            printf("%ld,",G->V[i].sibl[j]);
            //if (G->V[i].sibl[j] < 0 || G->V[i].sibl[j] >= G->n) printf("POZOR!!!"); getchar();}
        }
        printf("\n");*/
        if (G->V[i].edgecount > 1000)
        {
            printf("POZOR! Overflow!\n");
        }
    }
    //getchar();

    edges->clear();
    delete(edges);
}

// dealokujeme potrebnu pamat
void free_graph()
{
    unsigned long i;

    for (i=0;i<G->n;i++)
    {
        G->V[i].edgecount = 0;
        free(G->V[i].sibl);
    }

    free(G);
    G = NULL;

    //edges->clear();
}

bool are_adjacent(refer v, refer w)
{
    /*long i;
    for (i=0;i<G->V[v].edgecount;i++)
    {
        if (G->V[v].sibl[i] == w)
        {
            return true;
        }
    }*/
    return edges->isin(v,w);
}
