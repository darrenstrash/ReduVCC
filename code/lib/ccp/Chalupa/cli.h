#ifndef CLI_H
#define CLI_H
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <vector>
//#include <QTime>
#include "graphs.h"
#include "random_generator.h"
#include "algorithm_brelaz_ccp.h"
#include "algorithm_iggcc.h"
#include "algorithm_greedyindset.h"
#include "problem_ccp.h"
#include "time.h"
#include "statistics.h"
#include "timer.h"

class cli
{
private:
    long i,succ_rate;
////    unsigned long max_deg,min_deg;
    unsigned long long current_time,avg_time,avg_iter;
////    double avg_deg,stdev_deg;
    char filename[PATH_MAX];
    double avg_brelaz;
    double avg_iggcc;
    double avg_gis;
    double avg_rls_is;
    refer min_iggcc;
    unsigned long w,n_max;
    bool is_optimal;
//    QTime timer;
    FILE *source;
    graph G;
    refer max_label;
    int choose_instance(std::vector<std::vector<int>> const &adj_list, unsigned int num_v, unsigned long num_e);
    int generate_instance();
    int choose_algorithm(timer &t, double &time_to_solution, std::size_t clique_cover_offset);
    int compute_statistics();
    void sleep(unsigned long long milisec);
    void try_all_permutations();
    void generate_permutation(long permutation[], bool occupied[], refer level, algorithm *alg);
    refer* result_for_permutation; //[MAX_VERTICES];
    refer* histogram; //[MAX_VERTICES];
    refer* sizes; //[MAX_VERTICES];
    long   bad_scenario_counter;
////    long*  degree_distrib; //[MAX_VERTICES];
    problem_ccp *prob;
    double approx;
    long discrep;
    
    int seed;
    int mis;
    double t_elapsed;
    double t_limit;
    
public:
    std::vector<std::vector<int>> clique_cover;
    
    cli(int s, int m);
    ~cli();
    int start_cli(std::vector<std::vector<int>> const &adj_list, unsigned int num_v, unsigned long num_e, timer &total_timer, double &time_to_solution, double limit, std::size_t clique_cover_offset = 0);
};

#endif // CLI_H
