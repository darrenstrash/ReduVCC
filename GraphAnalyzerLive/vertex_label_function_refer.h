#ifndef VERTEX_LABEL_FUNCTION_REFER_H
#define VERTEX_LABEL_FUNCTION_REFER_H
#include "common.h"
#include "graphs.h"

class vertex_label_function_refer
{
private:
    refer vertices;
    refer k;
    refer (* vertex_label_matrix)[MAX_LABELS_CCP];
    class sorted_value_list
    {
    public:
        refer label;
        refer value;
    };
    class sorted_value_lists
    {
    public:
        refer sorted_values_count;
        sorted_value_list *sorted_values;
    };
    sorted_value_lists *lists;
public:
    vertex_label_function_refer(graph G, refer vertices, refer k);
    ~vertex_label_function_refer();
    void clear();
    refer get_value(refer vertex, refer label);
    void set_value(refer vertex, refer label, refer value);
};

#endif // VERTEX_LABEL_FUNCTION_SHORT_H
