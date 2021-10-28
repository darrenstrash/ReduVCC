
/******************************************************************************
 * vertex_queue.h
 * * Vertex queue for cascading reductions
 *
 *****************************************************************************/
 #ifndef QUEUE
 #define QUEUE


#include "data_structure/graph_access.h"
#include "redu_vcc.h"

class vertex_queue {
private:
    std::vector<bool> queue_status;
    std::vector<NodeID> queue;
public:

    vertex_queue(redu_vcc &reduVCC) {
      queue_status.assign(reduVCC.num_nodes, false);
    };
    ~vertex_queue() {};

    void init(redu_vcc &reduVCC) {
      for (NodeID v; v < reduVCC.num_nodes; v++) {
        queue_status[v] = true;
        queue.push_back(v);
      }
    };

    unsigned int size() { return queue.size(); };

    bool empty() { return queue.size() == 0; };
    void push(NodeID v) {
      if (!queue_status[v]) {
        queue.push_back(v);
        queue_status[v] = true;
      }
    };
    NodeID pop() {
      NodeID v = queue.back();
      queue.pop_back();
      queue_status[v] = false;
      return v;
    };

    void adjust_queue(redu_vcc &reduVCC, NodeID v){

      std::vector<std::vector<NodeID>> &adj_list = reduVCC.adj_list;
      std::vector<bool> &node_status = reduVCC.node_status;

      for (NodeID u : adj_list[v]) {
        if (!node_status[u]) continue;
        push(u);
      }
    };
};

#endif
