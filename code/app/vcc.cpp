
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

std::vector<bool> scratch1;
std::vector<bool> scratch2;

class Reduction {
public:
    static void removeVertex(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v);

    static void addCliqueToCover(std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<NodeID> &clique);
    static bool isSubset(std::vector<NodeID> &clique, NodeID &node, std::vector<NodeID> &neighborhood);


    virtual void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node) = 0;
    virtual void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) = 0;

};


void Reduction::removeVertex(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v){

    node_status[v] = false;
    for (NodeID u : adj_list[v]){
        for (unsigned int i = 0; i < adj_list[u].size(); i++){
            NodeID z = adj_list[u][i];
            if (z == v){
                adj_list[u].erase(adj_list[u].begin() + i);
                continue;
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

bool Reduction::isSubset(std::vector<NodeID> &clique, NodeID &node, std::vector<NodeID> &neighborhood){

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


class ISOReduction: public Reduction {
private:

	NodeID v;
    void addISOClique(std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
    void removeISO(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status);
    void removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a);

public:
    static bool validISO(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v);
    static bool validNeighbor(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {return;}
};

bool ISOReduction::validNeighbor(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u){
    // checks if u shares an edge with all w \in N(v)

    unsigned int i = 0;
    unsigned int j = 0;

    while (i < adj_list[v].size()) {
        NodeID x = adj_list[v][i];
        NodeID y = adj_list[u][j];

        if (x == u) {i++; continue;}

        if (j == adj_list[u].size()){return false;}

        if (x == y){i++; j++;}
        else if (x > y) {j++;}
        else {return false;}
    }

    return true;
}

bool ISOReduction::validISO(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v){
    // checks if v is an isolated vertex

    for (NodeID u : adj_list[v]) {
        bool valid_neighbor = ISOReduction::validNeighbor(adj_list, v, u);
        if (!valid_neighbor){return false;}
    }
    return true;

}

void ISOReduction::removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a){
    
    node_status[a] = false;
    for (NodeID b : adj_list[a]){
        for (unsigned int i = 0; i < adj_list[b].size(); i++){
            NodeID c = adj_list[b][i];
            if (c == a){
                adj_list[b].erase(adj_list[b].begin() + i);
                continue;
            }
        }
    }
    adj_list[a].erase(adj_list[a].begin(), adj_list[a].end());
}

void ISOReduction::removeISO(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status) {

    while (adj_list[v].size() > 0) {
        NodeID u = adj_list[v][0];
        // removeVertex(adj_list, node_status, q, u);
        removeNode(adj_list, node_status, u);
    }
    // removeVertex(adj_list, node_status, q, v);
    removeNode(adj_list, node_status, v);
}

void ISOReduction::addISOClique(std::vector<std::vector<NodeID>> &adj_list, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){

    std::vector<NodeID> clique;
    clique.push_back(v);

    for (NodeID u : adj_list[v]){
        clique.push_back(u);
    }

    addCliqueToCover(clique_cover, node_clique, clique);
}

void ISOReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node){

	v = node;
    // std::cout << "Reducing ISO... " << std::endl;
    addISOClique(adj_list, clique_cover, node_clique);
    removeISO(adj_list, node_status);
}


class D2Reduction: public Reduction {
private:
    NodeID v;
    NodeID u;
    NodeID w;

    std::vector<NodeID> N_u; // original neighborhood of u (excluding v)

    void removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a);
    void foldD2(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status);

public:
    static void assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w);
    static bool isTriangle(std::vector<std::vector<NodeID>> &adj_list, NodeID &u, NodeID &w);
    static bool validNeighbors(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID & w);
    static bool validD2(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};

void D2Reduction::assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w){

    NodeID x = adj_list[v][0];
    NodeID y = adj_list[v][1];

    if (adj_list[x].size() >= adj_list[y].size()) {u = y; w = x;}
    else {u = x; w = y;}
}

bool D2Reduction::isTriangle(std::vector<std::vector<NodeID>> &adj_list, NodeID &u, NodeID &w){

    for (NodeID x : adj_list[u]){
        if (x == w) {return true;}
    }
    return false;
}

bool D2Reduction::validNeighbors(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID & w){

    // scratch1 = N[u] / N[w]
    // scratch 2 = N[w] / N[u]
    for (NodeID x : adj_list[u]) {scratch1[x] = true;}
    for (NodeID x : adj_list[w]) {scratch2[x] = true; scratch1[x] = false;}
    for (NodeID x : adj_list[u]) {scratch2[x] = false;}

    for (NodeID x : adj_list[u]){
        if (!scratch1[x]) {continue;}
        for (NodeID y : adj_list[x]) {
            if (scratch2[y]) {return false;}
        }
    }

    return true;
}

bool D2Reduction::validD2(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v){

    if (adj_list[v].size() != 2) {return false;}

    NodeID u; NodeID w;
    D2Reduction::assignNodes(adj_list, v, u, w);

    if (D2Reduction::isTriangle(adj_list, u, w)) {return false;}

    if (D2Reduction::validNeighbors(adj_list, v, u, w)) {
        for (NodeID x : adj_list[u]) {scratch1[x] = false;}
        for (NodeID x : adj_list[w]) {scratch2[x] = false;}
        return true;} 

    for (NodeID x : adj_list[u]) {scratch1[x] = false;}
    for (NodeID x : adj_list[w]) {scratch2[x] = false;}
    return false;

}

void D2Reduction::removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a){

    node_status[a] = false;
    for (NodeID b : adj_list[a]){
        for (unsigned int i = 0; i < adj_list[b].size(); i++){
            NodeID c = adj_list[b][i];
            if (c == a){
                adj_list[b].erase(adj_list[b].begin() + i);
                continue;
            }
        }
    }
    adj_list[a].erase(adj_list[a].begin(), adj_list[a].end());
}

void D2Reduction::foldD2(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status) {

    for (NodeID x : adj_list[w]) {scratch1[x] = true;}

    for (NodeID x : adj_list[u]) {
        if (scratch1[x]) {continue;}
        adj_list[w].push_back(x);
        adj_list[x].push_back(w);
        std::sort(adj_list[x].begin(), adj_list[x].end());
    }
    std::sort(adj_list[w].begin(), adj_list[w].end());

    for (NodeID x: adj_list[w]) {
        scratch1[x] = false;
    }

    removeNode(adj_list, node_status, u);

}

void D2Reduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node){

    v = node;
    D2Reduction::assignNodes(adj_list, v, u, w);

    removeNode(adj_list, node_status, v);

    N_u = adj_list[u];
    foldD2(adj_list, node_status);
}


void D2Reduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){
    // std::cout << "Unreducing D2... " << std::endl;

    unsigned int fold_cliqueID = node_clique[u];
    std::vector<NodeID> fold_clique = clique_cover[fold_cliqueID];

    NodeID x = u;   // x denotes vertex connected to external
    NodeID y = w;   // y is connected to v

    for (NodeID p : N_u) {scratch1[p] = true;}
    for (NodeID p : fold_clique) {
        if (scratch1[p] == false) {
            x = w;
            y = u;
        }
    }
    for (NodeID p : N_u) {scratch1[p] = false;}

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


class TWINReduction: public Reduction {
private:
    NodeID v;
    NodeID u;
    NodeID w;
    NodeID x;
    NodeID y;

    std::vector<NodeID> N_w; // original neighborhood of u (excluding v)
    std::vector<NodeID> N_x;

    std::vector<NodeID> edge_nodes;
    NodeID nonedge_node;

    bool remove_type;

    bool removeType(std::vector<std::vector<NodeID>> &adj_list);

    void removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a);
    void removeTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
    void foldTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status);

public:
	static bool twinFound(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y);
    static void assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &w, NodeID &x, NodeID &y);
    static bool validTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u);

    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};

bool TWINReduction::twinFound(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &u, NodeID &w, NodeID &x, NodeID &y){

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

void TWINReduction::assignNodes(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &w, NodeID &x, NodeID &y){

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

bool TWINReduction::validTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u){

	if (adj_list[v].size() != 3) {return false;}

	NodeID w;
    NodeID x;
    NodeID y;

    TWINReduction::assignNodes(adj_list, v, w, x, y);

    if (TWINReduction::twinFound(adj_list, v, u, w, x, y)){
        return true;
    }

    return false;
}

void TWINReduction::removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a){

    node_status[a] = false;
    for (NodeID b : adj_list[a]){
        for (unsigned int i = 0; i < adj_list[b].size(); i++){
            NodeID c = adj_list[b][i];
            if (c == a){
                adj_list[b].erase(adj_list[b].begin() + i);
                continue;
            }
        }
    }
    adj_list[a].erase(adj_list[a].begin(), adj_list[a].end());
}


void TWINReduction::foldTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status){


	removeNode(adj_list, node_status, v);
    removeNode(adj_list, node_status, u);

    for (NodeID p : adj_list[y]){
        scratch1[p] = true;
    }

    for (NodeID p : adj_list[x]){

        if (p == y){
            continue;
        }
        if (scratch1[p]){
            continue;
        }
        adj_list[y].push_back(p);
        adj_list[p].push_back(y);

        std::sort(adj_list[p].begin(), adj_list[p].end());
    }

    for (NodeID p : adj_list[y]){
        scratch1[p] = true;
    }

    for (NodeID p : adj_list[w]){

        if (p == y){
            continue;
        }
        if (scratch1[p]){
            continue;
        }
        adj_list[y].push_back(p);
        adj_list[p].push_back(y);

        std::sort(adj_list[p].begin(), adj_list[p].end());
    }
    std::sort(adj_list[y].begin(), adj_list[y].end());

    for (NodeID z : adj_list[y]){
        scratch1[z] = false;
    }

    removeNode(adj_list, node_status, w);
    removeNode(adj_list, node_status, x);

 //    std::cout << "fold" << std::endl;
 //    std::cout << "y: " << y << std::endl;
	// std::cout << "N(y): ";
	// for (NodeID a : adj_list[y]){
	// 	std::cout << a << ", ";
	// }
	// std::cout << std::endl;
}


void TWINReduction::removeTWIN(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){

	std::vector<NodeID> clique1;
	clique1.push_back(v);
	removeNode(adj_list, node_status, v);

	for (NodeID p : edge_nodes){
        clique1.push_back(p);
    }
    std::sort(clique1.begin(), clique1.end());
    addCliqueToCover(clique_cover, node_clique, clique1);

    for (NodeID p : edge_nodes){
        removeNode(adj_list, node_status, p);
    }

    std::vector<NodeID> clique2;
    clique2.push_back(u);
    clique2.push_back(nonedge_node);
    std::sort(clique2.begin(), clique2.end());
    addCliqueToCover(clique_cover, node_clique, clique2);

    for (NodeID p : clique2){
        removeNode(adj_list, node_status, p);
    }

}

bool TWINReduction::removeType(std::vector<std::vector<NodeID>> &adj_list){

	for (NodeID p : adj_list[w]){
        if (p == x){
            edge_nodes.push_back(w);
            edge_nodes.push_back(x);
            nonedge_node = y;
            return true;
        }
        if (p == y){
            edge_nodes.push_back(w);
            edge_nodes.push_back(y);
            nonedge_node = x;
            return true;
        }
    }
    for (NodeID p : adj_list[x]){
        if (p == y){
            edge_nodes.push_back(x);
            edge_nodes.push_back(y);
            nonedge_node = w;
            return true;
        }
    }
    return false;
}

void TWINReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node){

	v = node;
	u = other_node;

	TWINReduction::assignNodes(adj_list, v, w, x, y);

	if (removeType(adj_list)){
		// std::cout << "Remove" << std::endl;
		remove_type = true;
		removeTWIN(adj_list, node_status, clique_cover, node_clique);
	}
	else {
		remove_type = false;
		N_w = adj_list[w];
		N_x = adj_list[x];
		foldTWIN(adj_list, node_status);
	}
}


void TWINReduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique){
	if (remove_type) {return;}

	NodeID p = y;
    NodeID q = x;
    NodeID r = w;

    unsigned int fold_cliqueID = node_clique[y];
    std::vector<NodeID> fold_clique = clique_cover[fold_cliqueID];

	if (isSubset(fold_clique, y, N_w)){
        p = w;
        q = x;
        r = y;
    }

    else if (isSubset(fold_clique, y, N_x)){
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


class DOMReduction: public Reduction {
private:

	NodeID v;
	NodeID u;

    void removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a);

public:
	static bool validDOM(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u);
    static bool nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &p);
    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique);
};

bool DOMReduction::nodeDominates(std::vector<std::vector<NodeID>> &adj_list, NodeID &v, NodeID &p){

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

bool DOMReduction::validDOM(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &v, NodeID &u){

	if (adj_list[v].size() > 2) {
        for (NodeID p : adj_list[v]){
            if (adj_list[v].size() < adj_list[p].size()){
                continue;
            }
            if (DOMReduction::nodeDominates(adj_list, v, p)){
                u = p;
                return true;
            }
        }
    }
    return false;
}

void DOMReduction::removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a){
    
    node_status[a] = false;
    for (NodeID b : adj_list[a]){
        for (unsigned int i = 0; i < adj_list[b].size(); i++){
            NodeID c = adj_list[b][i];
            if (c == a){
                adj_list[b].erase(adj_list[b].begin() + i);
                continue;
            }
        }
    }
    adj_list[a].erase(adj_list[a].begin(), adj_list[a].end());
}

void DOMReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &new_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node){

    v = node;
    u = other_node;

    removeNode(adj_list, node_status, v);
}

void DOMReduction::unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {

    unsigned int clique_ID = node_clique[u];
    std::vector<NodeID> clique = clique_cover[clique_ID];

    clique.push_back(v);
    std::sort(clique.begin(), clique.end());
    node_clique[v] = clique_ID;

    clique_cover[clique_ID] = clique;
}

class CROWNReduction: public Reduction {
private:

    void removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a);
	void addCrownCliques (std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &crown_cliques);

public:
    void reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node);
    void unreduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique) {return;}
};

void CROWNReduction::removeNode(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, NodeID &a){
    
    node_status[a] = false;
    for (NodeID b : adj_list[a]){
        for (unsigned int i = 0; i < adj_list[b].size(); i++){
            NodeID c = adj_list[b][i];
            if (c == a){
                adj_list[b].erase(adj_list[b].begin() + i);
                continue;
            }
        }
    }
    adj_list[a].erase(adj_list[a].begin(), adj_list[a].end());
}

void CROWNReduction::addCrownCliques (std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, std::vector<std::vector<int>> &crown_cliques) {
    
    // std::cout << crown_cliques.size() << std::endl;
    for (unsigned int i = 0; i < crown_cliques.size(); i++){
        std::vector<NodeID> clique;
        
        for (unsigned int j = 0; j < crown_cliques[i].size(); j++){
            int v = crown_cliques[i][j];
            NodeID old_v = new_to_old_map[v];
            
            removeNode(adj_list, node_status, old_v);
            clique.push_back(old_v);
        }
        std::sort(clique.begin(), clique.end());
        
        addCliqueToCover(clique_cover, node_clique, clique);
    }
}

void CROWNReduction::reduce(std::vector<std::vector<NodeID>> &adj_list, std::vector<bool> &node_status, std::vector<std::vector<int>> &int_adj_list, std::vector<NodeID> &new_to_old_map, std::vector<std::vector<NodeID>> &clique_cover, std::vector<unsigned int> &node_clique, NodeID &node, NodeID &other_node){

	branch_and_reduce_algorithm b_and_r(int_adj_list, int_adj_list.size());
	// std::cout << "searched for crown" << std::endl;
    if (b_and_r.lpCrownReduction()){
    	// std::cout << "crown found" << std::endl;
        addCrownCliques(adj_list, node_status, int_adj_list, new_to_old_map, clique_cover, node_clique, b_and_r.crown_cliques);
    }
}


class Reducer {
private:
    std::vector<Reduction*> reduction_stack;

    std::vector<std::vector<int>> new_adj_list;
    std::vector<int> old_to_new_map;
    std::vector<NodeID> new_to_old_map;

    void generateAdjList(graph_access &G);
    void generateNewAdjList(graph_access &G);
    unsigned int assignMaps(graph_access &G); 
    void addKernelCliques(std::vector<std::vector<int>> &clique_set);

    void performIsolatedReductions(graph_access &G);
    void performDegreeTwoReductions(graph_access &G);
    void performTwinReductions(graph_access &G);
    void performDomReductions(graph_access &G);
    void performCrownReductions(graph_access &G);


public:
    std::vector<std::vector<NodeID>> adj_list;
    std::vector<bool> node_status;

    std::vector<std::vector<NodeID>> clique_cover;
    std::vector<unsigned int> node_clique;

    Reducer(graph_access &G);
    void performReductions(graph_access &G);
    void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);
    void unwindReductions(graph_access &G);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
    unsigned int kernelSize(graph_access &G);

};

void Reducer::generateAdjList(graph_access &G){
    /* Generates adjacency list for graph */

    forall_nodes(G, v){
        std::vector<NodeID> N_v;    // open neighborhood of v

        forall_out_edges(G, e, v){
            NodeID u = G.getEdgeTarget(e);
            N_v.push_back(u);

        } endfor
        adj_list.push_back(N_v);

    } endfor
}

Reducer::Reducer(graph_access &G){

    // std::cout << "Creating reducer... " << std::endl;
    
    // reform graph
    generateAdjList(G);
    node_status.assign(G.number_of_nodes(), true);

    // set up clique cover
    node_clique.resize(G.number_of_nodes());
}

unsigned int Reducer::kernelSize(graph_access &G){

    unsigned int kernel_size = 0;

    forall_nodes(G, v) {
        if (node_status[v]){
            kernel_size++;
        }
    } endfor

    return kernel_size;
}


void Reducer::performIsolatedReductions(graph_access &G){
   
     bool vertexReduced = true;
    
     while (vertexReduced){
         vertexReduced = false;
         
         forall_nodes(G, v){
             NodeID u;
             
             if (!node_status[v]) {continue;}
             
             if (ISOReduction::validISO(adj_list, node_status, v)){
                vertexReduced = true;
                 Reduction *pReduction = nullptr;
                 pReduction = new ISOReduction();
                 pReduction->reduce(adj_list, node_status, new_adj_list, new_to_old_map, clique_cover, node_clique, v, u);
            }
         } endfor
    }
}

void Reducer::performDegreeTwoReductions(graph_access &G){


    bool vertexReduced = true;
    while (vertexReduced){
    vertexReduced = false; 
    forall_nodes(G, v){
    	NodeID u;
        if (!node_status[v]) {continue;}
        if (D2Reduction::validD2(adj_list, node_status, v)){
            vertexReduced = true;
            Reduction *pReduction = nullptr;
            pReduction = new D2Reduction();
            pReduction->reduce(adj_list, node_status, new_adj_list, new_to_old_map, clique_cover, node_clique, v, u);
        }
    } endfor
    }

}

void Reducer::performTwinReductions(graph_access &G){


    bool vertexReduced = true;
    while (vertexReduced){
    vertexReduced = false; 
    forall_nodes(G, v){
    	NodeID u;
        if (!node_status[v]) {continue;}
        if (TWINReduction::validTWIN(adj_list, node_status, v, u)){
            // std::cout << v << std::endl;
            vertexReduced = true;
            Reduction *pReduction = nullptr;
            pReduction = new TWINReduction();
            pReduction->reduce(adj_list, node_status,new_adj_list, new_to_old_map, clique_cover, node_clique, v, u);
        }
    } endfor
    }

}


void Reducer::performDomReductions(graph_access &G){

    bool vertexReduced = true;
    while (vertexReduced){
    vertexReduced = false; 
    forall_nodes(G, v){
    	NodeID u;
        if (!node_status[v]) {continue;}
        if (DOMReduction::validDOM(adj_list, node_status, v, u)){
            // if (v == 72430){
            //     std::cout << "REDUCED" << std::endl;
            // }
            vertexReduced = true;
            Reduction *pReduction = nullptr;
            pReduction = new DOMReduction();
            pReduction->reduce(adj_list, node_status, new_adj_list, new_to_old_map, clique_cover, node_clique, v, u);
        }
    } endfor
    }
}

void Reducer::performCrownReductions(graph_access &G){

    NodeID v; NodeID u;

    Reduction *pReduction = nullptr;
    pReduction = new CROWNReduction();
    generateNewAdjList(G);
    pReduction->reduce(adj_list, node_status, new_adj_list, new_to_old_map, clique_cover, node_clique, v, u);
    reduction_stack.push_back(pReduction);

}


void Reducer::performReductions(graph_access &G){

    unsigned int old_size = G.number_of_nodes();
    unsigned int new_size = 0;

    while (new_size != old_size){
        old_size = new_size;

        performIsolatedReductions(G);
        performDegreeTwoReductions(G);
        performTwinReductions(G);
        performDomReductions(G);
        performCrownReductions(G);

        new_size = kernelSize(G);
    }
   // std::cout << "missed" << std::endl;
   // performDomReductions(G);
   //    performDomReductions(G);

}

unsigned int Reducer::assignMaps(graph_access &G) {

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
    
    return u;
}


void Reducer::generateNewAdjList(graph_access &G) {

    unsigned int kernel_size = assignMaps(G);

    new_adj_list.clear();
    new_adj_list.resize(kernel_size);
    
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
        new_adj_list[new_v] = adj;
    }

    // int edges = 0;
    // for (unsigned int i = 0; i < new_adj_list.size(); i++){
    //     for (NodeID x : new_adj_list[i]) {
    //         std::cout << x + 1 << " ";
    //         edges++;
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << edges << std::endl;
}

void Reducer::addKernelCliques(std::vector<std::vector<int>> &clique_set){

    for (unsigned int i = 0; i < clique_set.size(); i++){
        std::vector<NodeID> clique;
        
        for (unsigned int j = 0; j < clique_set[i].size(); j++){
            int v = clique_set[i][j];
            NodeID old_v = new_to_old_map[v];
            
            Reduction::removeVertex(adj_list, node_status, old_v);
            clique.push_back(old_v);
        }
        std::sort(clique.begin(), clique.end());
        
        Reduction::addCliqueToCover(clique_cover, node_clique, clique);
    }

}

void Reducer::solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t) {

    bool nodes_remaining = false;
    for (unsigned int i = 0; i < G.number_of_nodes(); i++){
        if (node_status[i]) {
            nodes_remaining = true;
            break;
        }
    }

    if (!nodes_remaining) {return;}

    generateNewAdjList(G);
    
    unsigned int num_nodes = 0;
    unsigned long num_edges = 0;
    
    for (unsigned int i = 0; i < new_adj_list.size(); i++){
        for (unsigned int j = 0; j < new_adj_list[i].size(); j++){
           num_edges++;
        }
        num_nodes++;
    }
    
    cli *cli_instance;
    cli_instance = new cli(partition_config.seed);
    cli_instance->start_cli(new_adj_list, num_nodes, num_edges, t.elapsed(), 60);

    if (cli_instance->clique_cover.size() != 0){
        addKernelCliques(cli_instance->clique_cover);
    }
    else {
        std::cout << "Chalupa's algorithm unable to solve in given time." << std::endl;
    }
    delete(cli_instance);


}

void Reducer::unwindReductions(graph_access &G) {

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

    scratch1.assign(G.number_of_nodes(), false);
    scratch2.assign(G.number_of_nodes(), false);

    Reducer R(G);
    R.performReductions(G);
    R.analyzeGraph(graph_filename, G, s);
    // R.solveKernel(G, partition_config, s);
    // R.unwindReductions(G);
    // R.analyzeGraph(graph_filename, G, s);

}
