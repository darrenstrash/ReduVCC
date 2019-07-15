#include "algorithm_iggcc.h"

algorithm_iggcc::algorithm_iggcc()
{
    /*for (unsigned long i=0; i<MAX_T_MAX; i++)
    {
        for (unsigned long j=0; j<5; j++)
        {
            t_cliques[i][j] = 0;
        }
        t_contributions[i] = 0;
    }*/
}

long algorithm_iggcc::greedy_indset(graph G, long permutation_indset[])
{
    int i,j;
    refer vertex,indset_size;

    for (i=0;i<G->n;i++)
    {
        have_adjacent_indset_member[i] = false;
    }

    indset_size = 0;

    for (i=0;i<G->n;i++)
    {
        vertex = permutation_indset[i];
        // now we look, whether the next guy is not adjacent to some of the indset members
        if (! have_adjacent_indset_member[vertex])
        {
            indset_size++;
            for (j=0;j<G->V[vertex].edgecount;j++)
            {
                have_adjacent_indset_member[G->V[vertex].sibl[j]] = true;
            }
        }
    }

    return indset_size;
}

bool algorithm_iggcc::iggcc_ccp(graph G, refer *result, refer *indset_size, refer *initial_indset, refer initial_indset_size)
{
    unsigned long long t,t_max,indset_it,t_stag,t_stag_max;
    long remainder;
    long q,r,s,aux,w;
    long colors_count_old,colors_count,rand_handle;

    my_t_data = get_t_data();

    this->G = G;
    this->result = result;
    colors = G->n;
    t_max = 5000000;
    t_stag_max = 5000000;
    // t_stag_max = 1000000000;
    indset_it = 5;

    // initial permutation generation
    for (q=0;q<G->n;q++)
    {
        possibilities[q] = q;
    }
    for (q=0;q<G->n;q++)
    {
        r = generator.random(0,(G->n-1)-q);
        permutation[q] = possibilities[r];
        possibilities[r] = possibilities[(G->n-1)-q];
    }
    /*permutation[0] = 8; permutation[1] = 9;
    permutation[2] = 1; permutation[3] = 6;
    permutation[4] = 0; permutation[5] = 5;
    permutation[6] = 3; permutation[7] = 4;
    permutation[8] = 2; permutation[9] = 7;*/

    colors_count = generate_state_inverse_greedy(permutation);
    // indset generation
    if (NULL == initial_indset)
    {
        for (q=0;q<G->n;q++)
        {
            possibilities[q] = q;
        }
        for (q=0;q<G->n;q++)
        {
            r = generator.random(0,(G->n-1)-q);
            permutation_indset[q] = possibilities[r];
            possibilities[r] = possibilities[(G->n-1)-q];
        }
        fitness_indset = greedy_indset(G,permutation_indset);
    }
    else
    {
        for (q=0;q<G->n;q++)
        {
            free_indset_vertices[q] = true;
        }
        for (q=0;q<initial_indset_size;q++)
        {
            permutation_indset[q] = initial_indset[q];
            free_indset_vertices[initial_indset[q]] = false;
        }
        r = 0;
        for (q=0;q<G->n;q++)
        {
            if (free_indset_vertices[q])
            {
                possibilities[r] = q;
                r++;
            }
        }
        for (q=0;q<G->n-initial_indset_size;q++)
        {
            r = generator.random(0,(G->n-initial_indset_size-1)-q);
            permutation_indset[q+initial_indset_size] = possibilities[r];
            possibilities[r] = possibilities[(G->n-initial_indset_size-1)-q];
        }
        fitness_indset = greedy_indset(G,permutation_indset);
    }

    // evolution of the permutation
    t = 0;
    t_stag = 0;
    // while (t < t_max && fitness_indset < colors_count)
    while (t_stag < t_stag_max && fitness_indset < colors_count)
    {
        colors_count_old = colors_count;

        // greedy clique covering
        colors_count = generate_state_inverse_greedy(permutation);
        for (q=0;q<G->n;q++)
        {
            permutation[q] = auxiliary_state[q];
        }

        // determining the sizes of the classes
        remainder = G->n;
        for (q=0;q<=colors_count-1;q++)
        {
            class_sizes[q] = class_beginnings[q+1] - class_beginnings[q];
            remainder -= class_sizes[q];
        }
        class_sizes[colors_count] = remainder;

        if (colors_count_old > colors_count)
        {
            t_stag = 0;
        }
        else
        {
            t_stag++;
        }

        //if (colors_count_old > colors_count)
        if (t % 20 == 0)
        {
            /*long loners = 0, edges = 0, triangles = 0, fourcliques = 0;
            for (long dz=1;dz<=colors_count;dz++)
            {
                if (class_sizes[dz] == 1)
                {
                    loners++;
                }                
                if (class_sizes[dz] == 2)
                {
                    edges++;                    
                }                
                if (class_sizes[dz] == 3)
                {                    
                    triangles++;                    
                }                
                if (class_sizes[dz] == 4)
                {
                    fourcliques++;
                }                
            }
            my_t_data->t_cliques[t][1]+=loners;
            my_t_data->t_cliques[t][2]+=edges;
            my_t_data->t_cliques[t][3]+=triangles;
            my_t_data->t_cliques[t][4]+=fourcliques;
            my_t_data->t_contributions[t]++;*/
            //printf("[%ld,%ld,%ld,%ld](%ld),",loners,edges,triangles,fourcliques,fitness_indset);
            printf("[%ld](%ld),",colors_count,fitness_indset);

        }

        // reordering the classes
        for (q=0;q<G->n;q++)
        {
            new_ordering[q] = 0;
        }

        rand_handle = generator.random(5,12);
        //rand_handle = 13;
        //rand_handle = 12;
        if (t > 0 && rand_handle < 5)
        {
            // initializing the vertex lists
            for (q=1;q<=colors_count;q++)
            {
                vertex_lists[q].occ = false;
                vertex_lists[q].next = NULL;
                vertex_lists[q].prev = NULL;
                vertex_lists[q].end = &vertex_lists[q];
            }
            max_size = -1;
            for (q=1;q<=colors_count;q++)
            {
                if (max_size < class_sizes[q])
                {
                    max_size = class_sizes[q];
                }
                // if it is the first in a list
                if (! vertex_lists[class_sizes[q]].occ)
                {
                    vertex_lists[class_sizes[q]].occ = true;
                    vertex_lists[class_sizes[q]].vertex = q;
                }
                else
                {
                    // adding to the end of the equal-size list
                    vertex_list *current_vertex;
                    current_vertex = vertex_lists[class_sizes[q]].end;
                    current_vertex->next = new vertex_list();
                    current_vertex->next->vertex = q;
                    current_vertex->next->next = NULL;
                    current_vertex->next->prev = current_vertex;
                    vertex_lists[class_sizes[q]].end = current_vertex->next;
                }
            }
            // finding the new ordering and deleting all the auxiliary stuff
            q = 1;
            for (r=max_size;r>=1;r--)
            {
                if (! vertex_lists[r].occ)
                {
                    continue;
                }
                vertex_list *current_vertex;
                current_vertex = vertex_lists[r].end;
                while (current_vertex->prev != NULL)
                {
                    new_ordering[q] = current_vertex->vertex;
                    q++;
                    current_vertex = current_vertex->prev;
                    delete(current_vertex->next);
                }
                new_ordering[q] = current_vertex->vertex;
                q++;
                current_vertex->occ = false;
                current_vertex->end = current_vertex;
            }
        }
        else if (rand_handle < 10)
        {
            // reverse ordering
            for (q=1;q<=colors_count;q++)
            {
                new_ordering[q] = colors_count-q+1;
            }
        }
        else if (rand_handle < 13)
        {
            // random shuffle
            for (q=1;q<=colors_count;q++)
            {
                possibilities[q] = q;
            }
            for (q=1;q<=colors_count;q++)
            {
                r = generator.random(1,colors_count+1-q);
                new_ordering[q] = possibilities[r];
                possibilities[r] = possibilities[colors_count+1-q];
            }
        }
        else
        {
            // random first, other shifted
            r = generator.random(1,colors_count);
            new_ordering[1] = r;
            for (q=2;q<=r;q++)
            {
                new_ordering[q] = q-1;
            }
            for (q=r+1;q<=colors_count;q++)
            {
                new_ordering[q] = q;
            }
        }

        // applying the new ordering
        s = 0;
        for (q=1;q<=colors_count;q++)
        {
            current_class = new_ordering[q];
            for (r=0;r<class_sizes[current_class];r++)
            {
                new_permutation[s+r] = permutation[class_beginnings[current_class]+r];
            }
            new_class_beginnings[current_class] = s;
            s += class_sizes[current_class];
        }

        // pure IG-GCC
        for (q=0;q<=G->n;q++)
        {
            permutation[q] = new_permutation[q];
        }

        // updating the independent set bound
        for (w=0;w<indset_it;w++)
        {
            r = generator.random(1,G->n-1);
            aux = permutation_indset[r];
            for (s=0;s<r;s++)
            {
                new_permutation_indset[s+1] = permutation_indset[s];
            }
            for (s=r+1;s<G->n;s++)
            {
                new_permutation_indset[s] = permutation_indset[s];
            }
            new_permutation_indset[0] = aux;

            new_fitness_indset = greedy_indset(G,new_permutation_indset);
            if (new_fitness_indset >= fitness_indset)
            {
                for (s=0;s<G->n;s++)
                {
                    permutation_indset[s] = new_permutation_indset[s];
                }
                if (fitness_indset < new_fitness_indset)
                {
                    fitness_indset = new_fitness_indset;
                }
            }
        }
        //if (t % 1 == 0) {printf("(%ld/%ld)",colors_count,fitness_indset);}
        // printing the number of triangles in the found solution

        t++;
    }

    printf(" [iterations: %llu][independent set: %ld]\n",t,fitness_indset);
    *indset_size = fitness_indset;

    /*unsigned long i;
    printf("1-cliques: "); putchar('(');
    for (i=0;i<t_max;i++)
    {
        if (my_t_data->t_contributions[i] == 0)
        {
            break;
        }
        printf("%0.1lf(%lu),", (double)my_t_data->t_cliques[i][1]/(double)my_t_data->t_contributions[i],my_t_data->t_contributions[i]);
    }
    putchar(')');putchar('\n');
    printf("2-cliques: "); putchar('(');
    for (i=0;i<t_max;i++)
    {
        if (my_t_data->t_contributions[i] == 0)
        {
            break;
        }
        printf("%0.1lf(%lu),", (double)my_t_data->t_cliques[i][2]/(double)my_t_data->t_contributions[i],my_t_data->t_contributions[i]);
    }
    putchar(')');putchar('\n');
    printf("3-cliques: "); putchar('(');
    for (i=0;i<t_max;i++)
    {
        if (my_t_data->t_contributions[i] == 0)
        {
            break;
        }
        printf("%0.1lf(%lu),", (double)my_t_data->t_cliques[i][3]/(double)my_t_data->t_contributions[i],my_t_data->t_contributions[i]);
    }
    putchar(')');putchar('\n');*/

    return (fitness_indset == colors_count);
}

bool algorithm_iggcc::gcc_ccp(graph G, refer *result, long *permutation)
{
    long colors_count;

    this->G = G;
    this->result = result;
    colors = G->n;

    colors_count = generate_state_inverse_greedy(permutation);

    return true;
}

unsigned long long algorithm_iggcc::get_iterations()
{
    return 0;
}

long algorithm_iggcc::generate_state_inverse_greedy(long permutation[])
{
    long i,j,l,q,vertex,color;
    long current_colors;

    for (j=0;j<G->n;j++)
    {
        clique_sizes[j] = 0;
        size_copies[j] = 0;
        changed_sizes[j] = false;
        result[j] = 0;
        internal_positions[j] = 0;
    }

    current_colors = 0;
    for (j=0;j<G->n;j++)
    {
        vertex = permutation[j];

        // we have the vertex, we find a color for it
        color = 0;
        for (l=0;l<G->V[vertex].edgecount;l++)
        {
            i = G->V[vertex].sibl[l];
            q = result[i];
            if (q > 0)
            {
                if (! changed_sizes[q])
                {
                    size_copies[q] = clique_sizes[q];
                    changed_sizes[q] = true;
                }
                clique_sizes[q]--;
                if (clique_sizes[q] == 0)
                {
                    if (color == 0 || color > q)
                    {
                        color = q;
                    }
                }
            }
        }
        if (color == 0 && current_colors < colors)
        {
            current_colors++;
            color = current_colors;
        }
        // restoration of the sizes
        for (l=0;l<G->V[vertex].edgecount;l++)
        {
            i = G->V[vertex].sibl[l];
            q = result[i];
            if (q > 0 && changed_sizes[q])
            {
                clique_sizes[q] = size_copies[q];
                changed_sizes[q] = false;
            }
        }
        // coloring the vertex
        result[vertex] = color;
        internal_positions[vertex] = clique_sizes[color];
        clique_sizes[color]++;
    }

    // we create the conflicts list
    conflicts = 0;
    for (i=0;i<G->n;i++)
    {
        if (result[i] == 0)
        {
            conflicting_vertices_list[conflicts] = i;
            positioning_list[i] = conflicts;
            conflicts++;
        }
    }

    transform_representation();

    return current_colors;
}

void algorithm_iggcc::transform_representation()
{
    long i,sum;

    for (i=0;i<=colors;i++)
    {
        class_sizes[i] = 0;
    }
    for (i=0;i<G->n;i++)
    {
        auxiliary_state[i] = 0;
    }

    // counting the colors and finding the beginnings
    for (i=0;i<G->n;i++)
    {
        class_sizes[result[i]]++;
    }
    sum = 0;
    for (i=0;i<=colors;i++)
    {
        class_beginnings[i] = sum;
        sum += class_sizes[i];
    }
    for (i=0;i<=colors;i++)
    {
        class_sizes[i] = 0;
    }

    // putting the values into the auxiliary state
    for (i=0;i<G->n;i++)
    {
        if (internal_positions[i] >= 0)
        {
            auxiliary_state[class_beginnings[result[i]]+internal_positions[i]] = i;
        }
        else
        {
            auxiliary_state[class_beginnings[result[i]]+class_sizes[result[i]]] = i;
        }
        class_sizes[result[i]]++;
    }
}

