#include "algorithm_brelaz_ccp.h"

algorithm_brelaz_ccp::algorithm_brelaz_ccp(graph G, refer vertices)
{

    priorities.resize(MAX_VERTICES, 0);
    priorities_secondary.resize(MAX_VERTICES, 0);

    //neighbor_color_matrix = new short[vertices][MAX_LABELS_CCP];
    vertex_label_function = new vertex_label_function_refer(G,G->n,MAX_LABELS_CCP);
}

algorithm_brelaz_ccp::~algorithm_brelaz_ccp()
{
    //delete(neighbor_color_matrix);
    vertex_label_function->clear();
    delete(vertex_label_function);
}

void algorithm_brelaz_ccp::brelaz_ccp(graph G, refer *result, refer max_labels)
{
    unsigned long i,j,c,l,vertex;
    long long priority;
    long long priority_secondary;
    long color;
    long current_colors;
    long selected_candidate;
    long candidates_count;

    for (i=0;i<G->n;i++)
    {
        priorities[i] = 0;
        priorities_secondary[i] = (long long)(G->V[i].edgecount);
        class_sizes[i] = 0;
    }

    /*for (vertex=0;vertex<G->n;vertex++)
    {
        for (c=0;c<=MAX_LABELS_CCP;c++)
        {
            neighbor_color_matrix[vertex][c] = 0;
        }
    }*/
    vertex_label_function->clear();

    candidates_count = 0;
    current_colors = 0;
    for (j=0;j<G->n;j++)
    {
        // we find out, how many candidates do we actually have
        priority = -1;
        priority_secondary = -1;
        for (i=0;i<G->n;i++)
        {
            if (priorities[i] > priority || (priorities[i] == priority && priorities_secondary[i] > priority_secondary))
            {
                candidates_count = 1;
                priority = priorities[i];
                priority_secondary = priorities_secondary[i];
            }
            else if (priorities[i] == priority && priorities_secondary[i] == priority_secondary)
            {
                candidates_count++;
            }
        }
        selected_candidate = generator.random(0,candidates_count-1);
        // and we pick one of those candidates
        for (i=0;i<G->n;i++)
        {
            if (priorities[i] == priority && priorities_secondary[i] == priority_secondary)
            {
                if (selected_candidate)
                {
                    selected_candidate--;
                }
                else
                {
                    vertex = i;
                    break;
                }
            }
        }        
        // we have the vertex, we find a color for it
        color = 0;
        for (c=1;c<=max_labels;c++)
        {
            //if ((class_sizes[c] - neighbor_color_matrix[vertex][c]) == 0)
            if ((class_sizes[c] - (long)(vertex_label_function->get_value(vertex,c))) == 0)
            {
                color = c;
                break;
            }
        }
        // labeling the vertex        
        priorities[vertex] = -1;
        result[vertex] = color;                
        // all the non-neighbors now do have a colored non-neighbor
        l = 0;
        for (i=0;i<G->n;i++)
        {
            if (i == vertex)
            {
                continue;
            }
            if (l >= G->V[vertex].edgecount || G->V[vertex].sibl[l] != i)
            {
                //if (priorities[i] >= 0 && (class_sizes[color] - neighbor_color_matrix[i][color]) == 0)
                if (priorities[i] >= 0 && (class_sizes[color] - vertex_label_function->get_value(i,color)) == 0)
                {
                    priorities[i]++;
                }
            }
            else if (l < G->V[vertex].edgecount)
            {
                l++;
            }
        }
        for (l=0;l<G->V[vertex].edgecount;l++)
        {
            //neighbor_color_matrix[G->V[vertex].sibl[l]][color]++;
            vertex_label_function->set_value(G->V[vertex].sibl[l],color,vertex_label_function->get_value(G->V[vertex].sibl[l],color)+1);
        }
        // we update the class sizes
        class_sizes[color]++;
    }
}
