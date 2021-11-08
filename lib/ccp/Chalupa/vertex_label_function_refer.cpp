#include "vertex_label_function_refer.h"

vertex_label_function_refer::vertex_label_function_refer(graph G, refer vertices, refer k)
{
    refer v;

    this->vertices = vertices;
    this->k = k;

    //vertex_label_matrix = new refer[vertices][MAX_LABELS_CCP];

    lists = new sorted_value_lists[G->n+1];
    for (v=0;v<G->n;v++)
    {
        lists[v].sorted_values = new sorted_value_list[G->V[v].edgecount];
        lists[v].sorted_values_count = 0;
    }

}

vertex_label_function_refer::~vertex_label_function_refer()
{
    refer v;

    //delete(vertex_label_matrix);

    for (v=0;v<vertices;v++)
    {
        delete(lists[v].sorted_values);
    }
    delete(lists);
}

void vertex_label_function_refer::clear()
{
    refer v;

    /*refer j;
    for (v=0;v<vertices;v++)
    {
        for (j=0;j<=k;j++)
        {
            vertex_label_matrix[v][j] = 0;
        }
    }*/

    for (v=0;v<vertices;v++)
    {
        lists[v].sorted_values_count = 0;
    }
}

refer vertex_label_function_refer::get_value(refer vertex, refer label)
{
    refer left,right,current;

    //return vertex_label_matrix[vertex][label];

    // if the list is empty, then it is a zero
    if (lists[vertex].sorted_values_count == 0)
    {
        return 0;
    }
    else
    {
        //for (current=0;current<lists[vertex].sorted_values_count;current++)
        //{
        //    if (lists[vertex].sorted_values[current].label == label)
        //    {
        //        return lists[vertex].sorted_values[current].value;
        //    }
        //}

        left = 0;
        right = lists[vertex].sorted_values_count-1;
        while (right >= left)
        {
            current = (left+right)/2;
            if (lists[vertex].sorted_values[current].label >= label)
            {
                if (lists[vertex].sorted_values[current].label == label)
                {
                    return lists[vertex].sorted_values[current].value;
                }
                if (current > 0)
                {
                    right = current-1;
                }
                else
                {
                    right = current;
                    break;
                }
            }
            else if (lists[vertex].sorted_values[current].label < label)
            {
                if (current < lists[vertex].sorted_values_count-1)
                {
                    left = current+1;
                }
                else
                {
                    left = current;
                    break;
                }
            }
        }

        return 0;
    }
}

void vertex_label_function_refer::set_value(refer vertex, refer label, refer value)
{
    refer i,left,right,current;

    //vertex_label_matrix[vertex][label] = value;

    // put the value into the right list so that it remains sorted
    if (lists[vertex].sorted_values_count == 0)
    {
        if (value > 0)
        {
            lists[vertex].sorted_values[0].label = label;
            lists[vertex].sorted_values[0].value = value;
            lists[vertex].sorted_values_count++;
        }
    }
    else
    {
        // search for the position
        //for (current=0;current<lists[vertex].sorted_values_count;current++)
        //{
        //    // previous must be < and current must >= to the label, which is set
        //    if (lists[vertex].sorted_values[current].label >= label &&
        //        (current == 0 || lists[vertex].sorted_values[current-1].label < label))
        //    {
        //       break;
        //    }
        //    if (current == lists[vertex].sorted_values_count-1)
        //    {
        //        break;
        //    }
        //}

        left = 0;
        right = lists[vertex].sorted_values_count-1;
        while (lists[vertex].sorted_values[current].label < label ||
               (current != 0 && lists[vertex].sorted_values[current-1].label >= label))
        {
            current = (left+right)/2;
            if (lists[vertex].sorted_values[current].label == label)
            {
                break;
            }
            else if (lists[vertex].sorted_values[current].label > label)
            {
                if (current > 0)
                {
                    right = current-1;
                }
                else
                {
                    right = current;
                    break;
                }
            }
            else if (lists[vertex].sorted_values[current].label < label)
            {
                if (current < lists[vertex].sorted_values_count-1)
                {
                    left = current+1;
                }
                else
                {
                    left = current;
                    break;
                }
            }
        }

        // case 1: if we are at the end and the last label is still smaller than the inserted one
        if (current == lists[vertex].sorted_values_count-1 && lists[vertex].sorted_values[current].label < label)
        {
            // we insert the value at the end
            if (value > 0)
            {
                lists[vertex].sorted_values[current+1].label = label;
                lists[vertex].sorted_values[current+1].value = value;
                lists[vertex].sorted_values_count++;
            }
        }
        // case 2: if we found the value
        else if (lists[vertex].sorted_values[current].label == label)
        {
            // case 2a: the new value is non-zero, then we update it
            if (value > 0)
            {
                lists[vertex].sorted_values[current].value = value;
            }
            // case 2b: the new value is zero, then we delete it
            else
            {
                for (i=current+1;i<lists[vertex].sorted_values_count;i++)
                {
                    lists[vertex].sorted_values[i-1].label = lists[vertex].sorted_values[i].label;
                    lists[vertex].sorted_values[i-1].value = lists[vertex].sorted_values[i].value;
                }
                lists[vertex].sorted_values_count--;
            }
        }
        // case 3: if we did not find the value, we insert it
        else if (lists[vertex].sorted_values[current].label > label)
        {
            if (value > 0)
            {
                // shifting
                for (i=lists[vertex].sorted_values_count;i>=current+1;i--)
                {
                    lists[vertex].sorted_values[i].label = lists[vertex].sorted_values[i-1].label;
                    lists[vertex].sorted_values[i].value = lists[vertex].sorted_values[i-1].value;
                }
                lists[vertex].sorted_values[current].label = label;
                lists[vertex].sorted_values[current].value = value;
                lists[vertex].sorted_values_count++;
            }
        }
    }
}
