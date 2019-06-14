/***************A***************************************************************
 * vcc.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
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

bool checkNeighbor(graph_access &G, std::vector<bool> &node_status, NodeID &v, NodeID &u) {
    
    EdgeID e_v = G.get_first_edge(v);
    EdgeID e_u = G.get_first_edge(u);
    
    bool not_found = false;
    
    while (e_v != G.get_first_invalid_edge(v)) {
        NodeID x = G.getEdgeTarget(e_v);
        NodeID w = G.getEdgeTarget(e_u);
        
        if (node_status[x] == false){
            e_v++;
            continue;
        }
        
        if (x == u) {
            e_v++;
            continue;
        }
        
        if (e_u == G.get_first_invalid_edge(u)){
            not_found = true;
            break;
        }
        
        if (x == w){
            e_v++;
            e_u++;
        }
        else if (x > w) {
            e_u++;
        }
        else {
            not_found = true;
            break;
        }
    }
    
    return not_found;
}

void getIsolated(graph_access &G, std::vector<bool> &node_status, std::vector<NodeID> &isolated_vertices) {

  // int num_isolated = 0;

  int isolated_found = true;

  while (isolated_found){

      int num_isolated = 0;
    
      forall_nodes(G, v) {
          bool isolated = true;

          if (node_status[v] == false){
            continue;
          }

          forall_out_edges(G, e, v){
              
            NodeID u = G.getEdgeTarget(e);
              
            if (node_status[u] == false) {
                continue;
            }

              bool not_found = checkNeighbor(G, node_status, v, u);

            if (not_found){
                isolated = false;
                break;
            }
          } endfor

          if (isolated) {
          //std::cout << v << std::endl;
              num_isolated++;
              node_status[v] = false;
              
              forall_out_edges(G, e, v){
                  NodeID u = G.getEdgeTarget(e);
                  node_status[u] = false;
              } endfor
          
              isolated_vertices.push_back(v);
          }
        } endfor

        std::cout << num_isolated << std::endl;

        if (num_isolated == 0){
          isolated_found = false;
        }
    }
  
 

      //  std::cout << std::endl;
      
      //  return num_isolated;
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

  unsigned int num_nodes = G.number_of_nodes();
  
  std::vector<bool> node_status(num_nodes, true);

  //forall_nodes(G, v) {
  //vertex_status.push_back(true);
  //} endfor

  std::vector<NodeID> isolated_vertices;

  getIsolated(G, node_status, isolated_vertices);
    
    int vertices_remaining = 0;
    
    for (unsigned int i = 0; i < num_nodes; i++) {
        if (node_status[i] == true) {
            vertices_remaining++;
        }
    }

  for (unsigned int i = 0; i < isolated_vertices.size(); i++) {
      std::cout<< isolated_vertices[i] << ' ';
  }
    
    std::cout << std::endl << "Total nodes: " << num_nodes << std::endl;

    std::cout << "Number of cliques removed: " << isolated_vertices.size() << std::endl;
    
    std::cout << "Vertices remaining: " << vertices_remaining << std::endl;
  
  //  int num_isolated = getNumIsolated(G);
  //	std::cout << "Number of isolated vertices in G: " << num_isolated << std::endl;
}

