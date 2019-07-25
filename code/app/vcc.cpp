
/******************************************************************************
 * vcc.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "data_structure/matrix/online_distance_matrix.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

#include "mis/kernel/branch_and_reduce_algorithm.h"

#include "ccp/Chalupa/cli.h"
#include "ccp/Chalupa/random_generator.h"
#include <time.h>

class Queue {
private:
    unsigned int size;

    unsigned int start;
    unsigned int end;

    std::vector<NodeID> queue;
    std::vector<bool> node_in_queue;

public:
    Queue(graph_access &G, std::vector<bool> &node_status);
    bool empty() {return (start == end);}
    void push(NodeID &v);
    NodeID pop();
    bool in_queue(NodeID &v) {return node_in_queue[v];}
};

Queue::Queue(graph_access &G, std::vector<bool> &node_status){
    queue.resize(G.number_of_nodes() + 1);
    size = G.number_of_nodes() + 1;
    node_in_queue.assign(G.number_of_nodes(), false);

    start = 0;
    end = 0;

    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (!node_status[i]){
            continue;
        }
        queue[end] = i;
        node_in_queue[i] = true;
        end++;
    }

}

void Queue::push(NodeID &v) {
    end = end % size;
    queue[end] = v;
    end = (end + 1) % size;
    node_in_queue[v] = true;

}

NodeID Queue::pop(){
    NodeID v = queue[start];
    if (start == (size - 1)){
        end = end % size;
    }
    start = (start + 1) % size;
    node_in_queue[v] = false;
    return v;
}


class Reduction {
public:

    static void removeVertex(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v);
    static void addCliqueToCover(std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<NodeID> &clique);
    static void addIntCliquesToCover (std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &clique_set);

    bool cliqueIsSubset(std::vector<NodeID> &clique, NodeID &node, std::vector<NodeID> &neighborhood);

    virtual void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node) = 0;
    virtual void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) = 0;

};

void Reduction::removeVertex(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v){

    node_status[v] = false;
    for (NodeID u : adj_list[v]){
        for (unsigned int i = 0; i < adj_list[u].size(); i++){
            NodeID z = adj_list[u][i];
            if (z == v){
                adj_list[u].erase(adj_list[u].begin() + i);
                break;
            }
        }
    }
    adj_list[v].erase(adj_list[v].begin(), adj_list[v].end());
}

void Reduction::addCliqueToCover(std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<NodeID> &clique){

    std::sort(clique.begin(), clique.end());
    clique_cover.push_back(clique);
    unsigned int cliqueID = clique_cover.size() - 1;

    for (NodeID u : clique){
        node_clique[u] = cliqueID;
    }
}

bool Reduction::cliqueIsSubset(std::vector<NodeID> &clique, NodeID &node, std::vector<NodeID> &neighborhood){

    unsigned int i = 0;
    unsigned int j = 0;

    while (i < clique.size()){
        if (clique[i] == node){
            i++;
            continue;
        }
        if (j == neighborhood.size()){
            return false;
        }
        if (clique[i] < neighborhood[j]){
            return false;
        }
        if (clique[i] == neighborhood[j]){
            i++;
            j++;
        }
        else {
            j++;
        }
    }

    return true;
}


void Reduction::addIntCliquesToCover (std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &clique_set) {
    
    for (unsigned int i = 0; i < clique_set.size(); i++){
        std::vector<NodeID> clique;
        
        for (unsigned int j = 0; j < clique_set[i].size(); j++){
            int v = clique_set[i][j];
            NodeID old_v = new_to_old_map[v];
            
            removeVertex(adj_list, node_status, old_v);
            clique.push_back(old_v);
        }
        std::sort(clique.begin(), clique.end());
        
        addCliqueToCover(clique_cover, node_clique, clique);
    }
}

class IsolatedReduction: public Reduction {
private:
    void addIsolatedClique(std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &v);
    void addChangedNodesToQueue(std::vector<std::vector<NodeID>> &adj_list, Queue *q, NodeID &v);
    void removeIsolated(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, NodeID &v);

public:
    static bool isIsolated(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v);
    static bool neighborInIsoClique(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {return;}
};

bool IsolatedReduction::neighborInIsoClique(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u){

    unsigned int i = 0;
    unsigned int j = 0;

    while (i < adj_list[v].size()) {
        NodeID neighbor_v = adj_list[v][i];
        NodeID neighbor_u = adj_list[u][j];

        if (neighbor_v == u) {
            i++;
            continue;
        }

        if (j == adj_list[u].size()){
            return false;
        }

        if (neighbor_v == neighbor_u){
            i++;
            j++;
        }
        else if (neighbor_v > neighbor_u) {
            j++;
        }
        else {
            return false;
        }
    }

    return true;
}



bool IsolatedReduction::isIsolated(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v) {

    for (NodeID u : adj_list[v]) {
        bool valid_neighbor = IsolatedReduction::neighborInIsoClique(adj_list, v, u);
        if (!valid_neighbor){

            return false;
        }
    }
    return true;
}

void IsolatedReduction::addChangedNodesToQueue(std::vector<std::vector<NodeID>> &adj_list, Queue *q, NodeID &v){

    for (NodeID u : adj_list[v]){
        if (!q->in_queue(u)){
            q->push(u);
        }
    }
}

void IsolatedReduction::removeIsolated(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, NodeID &v) {

    while (adj_list[v].size() > 0) {
        NodeID u = adj_list[v][0];
        addChangedNodesToQueue(adj_list, q, u);
        removeVertex(adj_list, node_status, u);
    }
    removeVertex(adj_list, node_status, v);
}


void IsolatedReduction::addIsolatedClique(std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &v){

    std::vector<NodeID> clique;
    clique.push_back(v);

    for (NodeID u : adj_list[v]){
        clique.push_back(u);
    }

    addCliqueToCover(clique_cover, node_clique, clique);
}

void IsolatedReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &v, NodeID &partner_node){

//        addChangedNodesToQueue(adj_list, q, v);
    addIsolatedClique(adj_list, clique_cover, node_clique, v);
    removeIsolated(adj_list, node_status, q, v);
}


class DegreeTwoReduction: public Reduction {
private:
    NodeID v;
    NodeID u;
    NodeID w;

    std::vector<NodeID> original_neighborhood_u;

    void assignNeighborsBySize(std::vector<std::vector<NodeID>> &adj_list);

    void addChangedNodesToQueue(std::vector<std::vector<NodeID>> &adj_list, Queue *q);

    void foldDegreeTwo(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch);

public:
    static bool isDegreeTwo(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch, NodeID &v);
    static bool isTriangleClique(std::vector<std::vector<NodeID>> &adj_list, NodeID &v);
    static bool edgeBetweenNeighbors(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &scratch, NodeID &v, NodeID &u, NodeID &w);
//
    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};

bool DegreeTwoReduction::edgeBetweenNeighbors(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &scratch, NodeID &v, NodeID &u, NodeID &w){

//    for (NodeID x : adj_list[u]){
//        if (x == v || x == w){
//            continue;
//        }
//        for (NodeID y : adj_list[x]){
//            if (y == u || y == w){
//                continue;
//            }
//            if (scratch[y]) {
//                return true;
//            }
//        }
//    }
//    return false;
    
    for (unsigned int i = 0; i < adj_list[u].size(); i++){
        NodeID x = adj_list[u][i];
        
//        if (!node_status[x]){
//            continue;
//        }
        if (x == v || x == w) {
            continue;
        }
        for (unsigned int j = 0; j < adj_list[x].size(); j++){
            NodeID y = adj_list[x][j];
//            if (!node_status[y]){
//                continue;
//            }
            if (y == u || y == w){
                continue;
            }
            if (scratch[y]){
                return true;
            }
        }
    }
    return false;
}

bool DegreeTwoReduction::isTriangleClique(std::vector<std::vector<NodeID>> &adj_list, NodeID &v){

    NodeID u = adj_list[v][0];
    NodeID w = adj_list[v][1];

    if (adj_list[u].size() == 2 && adj_list[w].size() == 2){
        for (NodeID x : adj_list[u]){
            if (x == w){
                return true;
            }
        }
    }
    return false;
}


bool DegreeTwoReduction::isDegreeTwo(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch, NodeID &v){

    if (adj_list[v].size() != 2){
        return false;
    }

    if (isTriangleClique(adj_list, v)){
        return false;
    }

    NodeID u = adj_list[v][0];
    NodeID w = adj_list[v][1];

    for (NodeID x : adj_list[w]){
        scratch[x] = true;
    }

    for (NodeID x : adj_list[u]){
        scratch[x] = false;
    }

    bool not_valid = edgeBetweenNeighbors(adj_list, scratch, v, u, w);

    for (NodeID x : adj_list[w]){
        scratch[x] = false;
    }

    return !not_valid;
}


void DegreeTwoReduction::addChangedNodesToQueue (std::vector<std::vector<NodeID>> &adj_list, Queue *q) {


    if (!q->in_queue(w)){
        q->push(w);
    }
    for (NodeID p : adj_list[w]){
        if (!q->in_queue(p)){
            q->push(p);
        }
        for (NodeID x : adj_list[p]){
            if (!q->in_queue(x)){
                q->push(x);
            }
        }
    }

}

void DegreeTwoReduction::foldDegreeTwo(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch) {

    removeVertex(adj_list, node_status, v);
    
    for (NodeID z : adj_list[w]) {
        scratch[z] = true;
    }
    
    for (NodeID z : adj_list[u]) {
        if (z == w){
            continue;
        }
        if (scratch[z]){
            continue;
        }
        adj_list[w].push_back(z);
        adj_list[z].push_back(w);
        std::sort(adj_list[z].begin(), adj_list[z].end());
    }
    std::sort(adj_list[w].begin(), adj_list[w].end());
    
    for (NodeID z : adj_list[w]) {
        scratch[z] = false;
    }
    
    removeVertex(adj_list, node_status, u);
}


void DegreeTwoReduction::assignNeighborsBySize(std::vector<std::vector<NodeID>> &adj_list){

    NodeID x = adj_list[v][0];
    NodeID y = adj_list[v][1];

    if (adj_list[x].size() < adj_list[y].size()){
        u = x;
        w = y;

        return;
    }

    u = y;
    w = x;
}

void DegreeTwoReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node){

    v = node;

    assignNeighborsBySize(adj_list);

    original_neighborhood_u = adj_list[u];

    foldDegreeTwo(adj_list, node_status, scratch);
    
    addChangedNodesToQueue(adj_list, q);

}

void DegreeTwoReduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){

//    std::cout << "reduce d2" << std::endl;

    unsigned int fold_cliqueID = node_clique[u];
    std::vector<NodeID> fold_clique = clique_cover[fold_cliqueID];

    NodeID x = w;
    NodeID y = u;

    if (cliqueIsSubset(fold_clique, u, original_neighborhood_u)){
        x = u;
        y = w;
    }

    for (unsigned int i = 0; i < fold_clique.size(); i++){
        if (fold_clique[i] == u){
            fold_clique[i] = x;
            node_clique[x] = fold_cliqueID;
            std::sort(fold_clique.begin(), fold_clique.end());
            break;
        }
    }
    clique_cover[fold_cliqueID] = fold_clique;

    std::vector<NodeID> new_clique;
    new_clique.push_back(v);
    new_clique.push_back(y);
    std::sort(new_clique.begin(), new_clique.end());
    addCliqueToCover(clique_cover, node_clique, new_clique);
}


class TwinReduction: public Reduction {
private:
    NodeID v;
    NodeID u;
    NodeID w;
    NodeID x;
    NodeID y;

    bool remove_type_twin;

    std::vector<NodeID> original_neighborhood_w;
    std::vector<NodeID> original_neighborhood_x;

    std::vector<NodeID> edge_neighbors;
    NodeID non_edge_node;

    bool existsEdgeBetweenNeighbors (std::vector<std::vector<NodeID>> &adj_list);
    void removeTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
    void foldTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch);

    void addChangedNodesToQueueFoldType (std::vector<std::vector<NodeID>> &adj_list, Queue *q);
    void addChangedNodesToQueueRemoveType (std::vector<std::vector<NodeID>> &adj_list, Queue *q);

public:
    static bool isTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch, NodeID &v, NodeID &u);
    static void assignNeighborsBySize(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &w, NodeID &x, NodeID &y);
    static bool twinFound(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};

bool TwinReduction::twinFound(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y){

    for (NodeID p : adj_list[w]){
        if (p == v){
            continue;
        }
        if (adj_list[p].size() != 3){
            continue;
        }

        for (NodeID q : adj_list[x]){
            if (q > p){
                break;
            }
            if (q == p){
                for (NodeID r : adj_list[y]){
                    if (r > p){
                        break;
                    }
                    if (r == p){
                        u = p;
                        return true;
                    }
                }
                break;
            }
        }
    }
    return false;
}

void TwinReduction::assignNeighborsBySize(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &w, NodeID &x, NodeID &y){

    NodeID a = adj_list[v][0];
    NodeID b = adj_list[v][1];
    NodeID c = adj_list[v][2];

    if (adj_list[a].size() < adj_list[b].size() && adj_list[a].size() < adj_list[c].size()) {
        w = a;
        if (adj_list[b].size() < adj_list[c].size()){
            x = b;
            y = c;
        }
        else {
            x = c;
            y = b;
        }
    }
    else if (adj_list[b].size() < adj_list[a].size() && adj_list[b].size() < adj_list[c].size()) {
        w = b;
        if (adj_list[a].size() < adj_list[c].size()){
            x = a;
            y = c;
        }
        else {
            x = c;
            y = a;
        }
    }
    else {
        w = c;
        if (adj_list[a].size() < adj_list[b].size()){
            x = a;
            y = b;
        }
        else {
            x = b;
            y = a;
        }
    }
}

bool TwinReduction::isTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch, NodeID &v, NodeID &u){

    if (adj_list[v].size() != 3){
        return false;
    }

    NodeID w;
    NodeID x;
    NodeID y;

    TwinReduction::assignNeighborsBySize(adj_list, v, w, x, y);

    if (TwinReduction::twinFound(adj_list, v, u, w, x, y)){
        return true;
    }

    return false;
}

void TwinReduction::addChangedNodesToQueueFoldType (std::vector<std::vector<NodeID>> &adj_list, Queue *q) {


    if (!q->in_queue(y)){
        q->push(y);
    }
    for (NodeID p : adj_list[y]){
        if (!q->in_queue(p)){
            q->push(p);
        }
        for (NodeID r : adj_list[p]){
            if (!q->in_queue(r)){
                q->push(r);
            }
        }
    }

}

void TwinReduction::foldTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<bool> &scratch) {

    removeVertex(adj_list, node_status, v);
    removeVertex(adj_list, node_status, u);

    for (NodeID p : adj_list[y]){
        scratch[p] = true;
    }

    for (NodeID q : adj_list[x]){

        if (q == y){
            continue;
        }
        if (scratch[q]){
            continue;
        }
        adj_list[y].push_back(q);
        adj_list[q].push_back(y);

        std::sort(adj_list[q].begin(), adj_list[q].end());
    }

    for (NodeID r : adj_list[y]){
        scratch[r] = true;
    }

    for (NodeID p : adj_list[w]){

        if (p == y){
            continue;
        }
        if (scratch[p]){
            continue;
        }
        adj_list[y].push_back(p);
        adj_list[p].push_back(y);

        std::sort(adj_list[p].begin(), adj_list[p].end());
    }
    std::sort(adj_list[y].begin(), adj_list[y].end());


    for (NodeID z : adj_list[y]){
        scratch[z] = false;
    }

    removeVertex(adj_list, node_status, w);
    removeVertex(adj_list, node_status, x);

}


void TwinReduction::addChangedNodesToQueueRemoveType (std::vector<std::vector<NodeID>> &adj_list, Queue *q){


    for (NodeID p : adj_list[v]){

        if (!q->in_queue(p)){
            q->push(p);
        }
        for (NodeID r : adj_list[p]){
            if (!q->in_queue(r)){
                q->push(r);
            }
            for (NodeID s : adj_list[w]){
                if (!q->in_queue(s)){
                    q->push(s);
                }
            }
        }
    }
}

void TwinReduction::removeTwin(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){

    std::vector<NodeID> clique1;
    clique1.push_back(v);
    for (NodeID p : edge_neighbors){
        clique1.push_back(p);
    }
    std::sort(clique1.begin(), clique1.end());
    addCliqueToCover(clique_cover, node_clique, clique1);

    for (NodeID p : clique1){
        removeVertex(adj_list, node_status, p);
    }

    std::vector<NodeID> clique2;
    clique2.push_back(u);
    clique2.push_back(non_edge_node);
    std::sort(clique2.begin(), clique2.end());
    addCliqueToCover(clique_cover, node_clique, clique2);

    for (NodeID p : clique2){
        removeVertex(adj_list, node_status, p);
    }
}

bool TwinReduction::existsEdgeBetweenNeighbors(std::vector<std::vector<NodeID>> &adj_list){

    for (NodeID p : adj_list[w]){
        if (p == x){
            edge_neighbors.push_back(w);
            edge_neighbors.push_back(x);
            non_edge_node = y;
            return true;
        }
        if (p == y){
            edge_neighbors.push_back(w);
            edge_neighbors.push_back(y);
            non_edge_node = x;
            return true;
        }
    }
    for (NodeID p : adj_list[x]){
        if (p == y){
            edge_neighbors.push_back(x);
            edge_neighbors.push_back(y);
            non_edge_node = w;
            return true;
        }
    }
    return false;
}

void TwinReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node){

    v = node;
    u = partner_node;

    TwinReduction::assignNeighborsBySize(adj_list, v, w, x, y);

    remove_type_twin = existsEdgeBetweenNeighbors(adj_list);
    if (remove_type_twin){
        addChangedNodesToQueueRemoveType(adj_list, q);
        removeTwin(adj_list, node_status, clique_cover, node_clique);
    }
    else {
        original_neighborhood_w = adj_list[w];
        original_neighborhood_x = adj_list[x];
        foldTwin(adj_list, node_status, scratch);
        addChangedNodesToQueueFoldType(adj_list, q);
    }
}


void TwinReduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {

    if (remove_type_twin){
        return;
    }

    unsigned int fold_cliqueID = node_clique[y];
    std::vector<NodeID> fold_clique = clique_cover[fold_cliqueID];

    NodeID p = y;
    NodeID q = x;
    NodeID r = w;

    if (cliqueIsSubset(fold_clique, y, original_neighborhood_w)){
        p = w;
        q = x;
        r = y;
    }

    else if (cliqueIsSubset(fold_clique, y, original_neighborhood_x)){
        p = x;
        q = w;
        r = y;
    }

    for (unsigned int i = 0; i < fold_clique.size(); i++){
        if (fold_clique[i] == y){
            fold_clique[i] = p;
            node_clique[p] = fold_cliqueID;
            std::sort(fold_clique.begin(), fold_clique.end());
            break;
        }
    }
    clique_cover[fold_cliqueID] = fold_clique;

    std::vector<NodeID> new_clique1;
    new_clique1.push_back(v);
    new_clique1.push_back(q);
    std::sort(new_clique1.begin(), new_clique1.end());
    addCliqueToCover(clique_cover, node_clique, new_clique1);

    std::vector<NodeID> new_clique2;
    new_clique2.push_back(u);
    new_clique2.push_back(r);
    std::sort(new_clique2.begin(), new_clique2.end());
    addCliqueToCover(clique_cover, node_clique, new_clique2);
}


class DominationReduction: public Reduction {
private:
    NodeID v;
    NodeID u;

    void addChangedNodesToQueue(std::vector<std::vector<NodeID>> &G, Queue *q);


public:
    static bool isDominant(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u);
    static bool nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &p);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};


bool DominationReduction::nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &p){

    unsigned int i = 0;
    unsigned int j = 0;

    while (j < adj_list[p].size()){
        if (adj_list[p][j] == v){
            j++;
            continue;
        }
        if (i == adj_list[v].size()){
            return false;
        }
        if (adj_list[v][i] == p){
            i++;
            continue;
        }
        if (adj_list[v][i] > adj_list[p][j]){
            return false;
        }
        if (adj_list[v][i] < adj_list[p][j]){
            i++;
        }
        else {
            i++;
            j++;
        }
    }
    return true;
}


bool DominationReduction::isDominant(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u) {

    if (adj_list[v].size() > 2) {
        for (NodeID p : adj_list[v]){
            if (adj_list[v].size() < adj_list[p].size()){
                continue;
            }
            if (DominationReduction::nodeDominates(adj_list, v, p)){
                u = p;
                return true;
            }
        }
    }
    return false;
}


void DominationReduction::addChangedNodesToQueue(std::vector<std::vector<NodeID>> &adj_list, Queue *q){

    for (NodeID w : adj_list[v]){
        if (!q->in_queue(w)){
            q->push(w);
        }
    }
}


void DominationReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &node, NodeID &partner_node){

    v = node;
    u = partner_node;

    addChangedNodesToQueue(adj_list, q);
    removeVertex(adj_list, node_status, v);
}

void DominationReduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {

    unsigned int clique_ID = node_clique[u];
    std::vector<NodeID> clique = clique_cover[clique_ID];

    clique.push_back(v);
    std::sort(clique.begin(), clique.end());
    node_clique[v] = clique_ID;

    clique_cover[clique_ID] = clique;
}


class CrownReduction: public Reduction {
private:
    
    void addCrownCliquesToCover(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &crown_cliques);
    
public:
    
    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &v, NodeID &partner_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {return;}

    
};

void CrownReduction::addCrownCliquesToCover (std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &crown_cliques) {
    
    for (unsigned int i = 0; i < crown_cliques.size(); i++){
        std::vector<NodeID> clique;
        
        for (unsigned int j = 0; j < crown_cliques[i].size(); j++){
            int v = crown_cliques[i][j];
            NodeID old_v = new_to_old_map[v];
            
            removeVertex(adj_list, node_status, old_v);
            clique.push_back(old_v);
        }
        std::sort(clique.begin(), clique.end());
        
        addCliqueToCover(clique_cover, node_clique, clique);
    }
}

void CrownReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, Queue *q, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<bool> &scratch, NodeID &v, NodeID &partner_node) {
    
    branch_and_reduce_algorithm b_and_r(int_adj_list, int_adj_list.size());
    if (b_and_r.lpCrownReduction()){
        addCrownCliquesToCover(adj_list, node_status, int_adj_list, new_to_old_map, clique_cover, node_clique, b_and_r.crown_cliques);
    }
}


class Reducer {
private:
    Queue *q;
    std::vector<Reduction*> reduction_stack;

    std::vector<bool> bool_scratch;
    
    std::vector<std::vector<int>> int_adj_list;
    std::vector<int> old_to_new_map;
    std::vector<NodeID> new_to_old_map;

    void makeAdjList(graph_access &G);
    void remakeAdjList(graph_access &G);
    
    int assignAdjMaps(graph_access &G);

    void performIsolatedReductions(graph_access &G);
    void performDegreeTwoReductions(graph_access &G);
    void performTwinReductions(graph_access &G);
    void performDominationReductions(graph_access &G);
    void performCrownReductions(graph_access &G);

public:
    unsigned int graph_size;
    unsigned int remaining_nodes;

    std::vector<std::vector<NodeID>> adj_list;
    std::vector<bool> node_status;

    std::vector<std::vector<NodeID>> clique_cover;
    std::vector<unsigned int> node_clique;

    Reducer(graph_access &G);
    unsigned int getSubgraphSize(graph_access &G);
    void performReductions(graph_access &G);
    void computeSubgraph(graph_access &G, PartitionConfig &partition_config);
    void unwindReductions(graph_access &G);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
};

void Reducer::makeAdjList(graph_access &G){

    forall_nodes(G, v){
        std::vector<NodeID> neighbors_of_v;

        forall_out_edges(G, e, v){
            NodeID u = G.getEdgeTarget(e);
            neighbors_of_v.push_back(u);

        } endfor
        adj_list.push_back(neighbors_of_v);

    } endfor
}

int Reducer::assignAdjMaps(graph_access &G) {
    
    old_to_new_map.clear();
    old_to_new_map.resize(G.number_of_nodes());
    new_to_old_map.clear();
    new_to_old_map.resize(G.number_of_nodes());
    
    int u = 0;
    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (!node_status[i]){
            continue;
        }
        old_to_new_map[i] = u;
        new_to_old_map[u] = i;
        u++;
    }
    
    int sub_graph_size = u;
    return sub_graph_size;
}

void Reducer::remakeAdjList(graph_access &G) {
    
    int sub_graph_size = assignAdjMaps(G);
    
    int_adj_list.clear();
    int_adj_list.resize(sub_graph_size);
    
    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (!node_status[i]){
            continue;
        }
        int new_v = old_to_new_map[i];
        std::vector<int> adj;
        for (unsigned int j = 0; j < adj_list[i].size(); j++){
            NodeID u  = adj_list[i][j];
            int new_u = old_to_new_map[u];
            adj.push_back(new_u);
        }
        std::sort(adj.begin(), adj.end());
        int_adj_list[new_v] = adj;
    }
    
    
}

Reducer::Reducer(graph_access &G){

    graph_size = 0;
    remaining_nodes = G.number_of_nodes();

    makeAdjList(G);
    node_status.assign(G.number_of_nodes(), true);
    node_clique.resize(G.number_of_nodes());

    bool_scratch.assign(G.number_of_nodes(), false);
}


unsigned int Reducer::getSubgraphSize(graph_access &G){

    unsigned int nodes_in_graph = 0;

    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (node_status[i]){
            nodes_in_graph++;
        }
    }

    return nodes_in_graph;

}


void Reducer::performIsolatedReductions(graph_access &G){
    q = new Queue(G, node_status);
    for(;;){
        bool reductions_complete = q->empty();
        if (reductions_complete){
            break;
        }

        NodeID v = q->pop();

//        int vertices_reduced = 0;
//        forall_nodes(G, v){

            if (!node_status[v]){
                continue;
            }

            NodeID u;
            Reduction *pReduction = nullptr;

            if (IsolatedReduction::isIsolated(adj_list, node_status, v)){
//                vertices_reduced++;
                pReduction = new IsolatedReduction();
            }

            else {
                continue;
            }
            pReduction->reduce(adj_list, node_status, q, int_adj_list, new_to_old_map, clique_cover, node_clique, bool_scratch, v, u);
            reduction_stack.push_back(pReduction);

//        } endfor
//
//        if (vertices_reduced == 0){
//            break;
//        }
    }

    delete q;
}

void Reducer::performDegreeTwoReductions(graph_access &G){
    q = new Queue(G, node_status);
    for(;;){
        bool reductions_complete = q->empty();
        if (reductions_complete){
            break;
        }

        NodeID v = q->pop();

//                        int vertices_reduced = 0;
//                        forall_nodes(G, v){

        if (!node_status[v]){
            continue;
        }

        NodeID u;
        Reduction *pReduction = nullptr;

        if (DegreeTwoReduction::isDegreeTwo(adj_list, node_status, bool_scratch, v)){
//                                        vertices_reduced++;
            pReduction = new DegreeTwoReduction();
        }

        else {
            continue;
        }
        pReduction->reduce(adj_list, node_status, q, int_adj_list, new_to_old_map, clique_cover, node_clique, bool_scratch, v, u);
        reduction_stack.push_back(pReduction);

//                        } endfor
//                        if (vertices_reduced == 0){
//
//                            break;
//                        }
    }

    delete q;
}


void Reducer::performTwinReductions(graph_access &G){
    q = new Queue(G, node_status);
    for(;;){
                bool reductions_complete = q->empty();
                if (reductions_complete){
                    break;
                }

                NodeID v = q->pop();
        
//        int vertices_reduced = 0;
//        forall_nodes(G, v){
        
            if (!node_status[v]){
                continue;
            }
            
            NodeID u;
            Reduction *pReduction = nullptr;
            
            if (TwinReduction::isTwin(adj_list, node_status, bool_scratch, v, u)){
//                vertices_reduced++;
                pReduction = new TwinReduction();
            }
            
            else {
                continue;
            }
            pReduction->reduce(adj_list, node_status, q, int_adj_list, new_to_old_map, clique_cover, node_clique, bool_scratch, v, u);
            reduction_stack.push_back(pReduction);
            //
//        } endfor
//
//        if (vertices_reduced == 0){
//            break;
//        }
    }
    
    delete q;
}


void Reducer::performDominationReductions(graph_access &G){
    q = new Queue(G, node_status);
    for(;;){
                bool reductions_complete = q->empty();
                if (reductions_complete){
                    break;
                }

                NodeID v = q->pop();

//        int vertices_reduced = 0;
//        forall_nodes(G, v){

            if (!node_status[v]){
                continue;
            }

            NodeID u;
            Reduction *pReduction = nullptr;

            if (DominationReduction::isDominant(adj_list, node_status, v, u)){
//                vertices_reduced++;
                pReduction = new DominationReduction();
            }

            else {
                continue;
            }
            pReduction->reduce(adj_list, node_status, q, int_adj_list, new_to_old_map, clique_cover, node_clique, bool_scratch, v, u);
            reduction_stack.push_back(pReduction);

//        } endfor

//        if (vertices_reduced == 0){
//            break;
//        }
    }

    delete q;
}

void Reducer::performCrownReductions(graph_access &G){
    
    remakeAdjList(G);
    
    NodeID v;
    NodeID u;
    
    Reduction *pReduction = new  CrownReduction();
    pReduction->reduce(adj_list, node_status, q, int_adj_list, new_to_old_map, clique_cover, node_clique, bool_scratch, v, u);
    reduction_stack.push_back(pReduction);
}


void Reducer::performReductions(graph_access &G){

    unsigned int sub_graph_size = 0;

    while (sub_graph_size != remaining_nodes){
        sub_graph_size =  remaining_nodes;
        performIsolatedReductions(G);
        performDegreeTwoReductions(G);
        
        performCrownReductions(G);
        
        performTwinReductions(G);
        performDominationReductions(G);

        remaining_nodes = getSubgraphSize(G);
    }
}

void Reducer::computeSubgraph(graph_access &G, PartitionConfig &partition_config){
    
    if (remaining_nodes == 0){
        return;
    }
    
    remakeAdjList(G);
    
    unsigned int num_nodes = 0;
    unsigned long num_edges = 0;
    
    for (unsigned int i = 0; i < int_adj_list.size(); i++){
        for (unsigned int j = 0; j < int_adj_list[i].size(); j++){
            num_edges++;
        }
        num_nodes++;
    }
    
    cli *cli_instance;
    cli_instance = new cli(partition_config.seed);
    cli_instance->start_cli(int_adj_list, num_nodes, num_edges);

    
    Reduction::addIntCliquesToCover(adj_list, node_status, int_adj_list, new_to_old_map, clique_cover, node_clique, cli_instance->clique_cover);
    
    delete(cli_instance);
}


void Reducer::unwindReductions(graph_access &G){

    while (reduction_stack.size() > 0){
        Reduction *pReduction = reduction_stack.back();
        reduction_stack.pop_back();

        pReduction->unreduce(adj_list, node_status, clique_cover, node_clique);

        delete pReduction;
    }

}

void Reducer::analyzeGraph(std::string &filename, graph_access &G, timer &t){

    std::cout << filename << ", ";

    std::cout << G.number_of_nodes() << ", ";
    std::cout << G.number_of_edges() << ", ";

    unsigned int remaining_nodes = 0;
    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (node_status[i]) {
            remaining_nodes++;
        }
    }

    std::cout << remaining_nodes << ", ";

    std::cout << t.elapsed() << ", ";

    std::cout << clique_cover.size() << ", ";

    std::vector<bool> in_clique(G.number_of_nodes(), false);
    for (std::vector<NodeID> clique : clique_cover){
        for (NodeID v : clique){
            if(in_clique[v]){
                std::cout << "Duplicate: " << v << ", ";
            }
            else {
                in_clique[v] = true;
            }
        }
    }

    for (bool in : in_clique){
        if (!in){
            std::cout << "Error" << std::endl;
            return;
        }
    }

    std::cout << "Complete" << std::endl;
}


class CliqueCoverBranchAndReduce {
    
    void makeAdjMatrix(graph_access &G);
    
public:
    
    char** adjacency_matrix;
    
    CliqueCoverBranchAndReduce(graph_access &G);
};

void CliqueCoverBranchAndReduce::makeAdjMatrix(graph_access &G){
    
    
    
    unsigned int nodes = G.number_of_nodes();
    adjacency_matrix = (char **) malloc (sizeof(char **)*nodes);

    for (unsigned int i = 0; i < nodes; i++){
        adjacency_matrix[i] = (char *) malloc (sizeof(char)*nodes);
    }

    for (unsigned int i = 0; i < nodes; i++){
        for (unsigned int j = 0; j < nodes; j++){
            adjacency_matrix[i][j] = 0;
        }
    }

    forall_nodes(G, v){
        forall_out_edges(G, e, v){
            NodeID u = G.getEdgeTarget(e);
            adjacency_matrix[v][u] = 1;
        } endfor
    } endfor
}

CliqueCoverBranchAndReduce::CliqueCoverBranchAndReduce(graph_access &G){
    
    makeAdjMatrix(G);
}

int main(int argn, char **argv) {

    PartitionConfig partition_config;
    std::string graph_filename;

    bool is_graph_weighted = false;
    bool suppress_output   = false;
    bool recursive         = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if(ret_code) {
        return 0;
    }

    std::streambuf* backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if(suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }

    partition_config.LogDump(stdout);
    graph_access G;

    timer t;
    graph_io::readGraphWeighted(G, graph_filename);
    //  std::cout << "io time: " << t.elapsed()  << std::endl;

    timer s;
    std::vector<std::vector<NodeID>> adj_list;

    
    Reducer R(G);
//    R.performReductions(G);
    R.computeSubgraph(G, partition_config);
//    R.unwindReductions(G);
    R.analyzeGraph(graph_filename, G, s);
}
