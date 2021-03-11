
/******************************************************************************
 * redu_vcc.h
 * *
 *
 *****************************************************************************/

#ifndef REDU_VCC
#define REDU_VCC

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "ccp/Chalupa/cli.h"
#include <time.h>

#include "redu_structure.h"


class redu_vcc : public redu_structure {

  public:

    redu_vcc() {};
    redu_vcc(graph_access &G) : redu_structure(G) {};
    redu_vcc(graph_access &G, PartitionConfig &partition_config);
    virtual ~redu_vcc() {};

    void getMIS(std::string file);

    void analyzeGraph(std::string &filename, graph_access &G, timer &t);
    void solveKernel(graph_access &G, PartitionConfig &partition_config, timer &t);
    void writeKernel(graph_access &G, std::string &filename);


};

#endif
