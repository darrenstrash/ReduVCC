#include "cli.h"

cli::cli()
{
    result_for_permutation = new refer[MAX_VERTICES];
    histogram              = new refer[MAX_VERTICES];
    sizes                  = new refer[MAX_VERTICES];
    degree_distrib         = new long[MAX_VERTICES];
}

cli::~cli() 
{
    delete[] result_for_permutation;
    delete[] histogram;
    delete[] sizes;
    delete[] degree_distrib;
}

void cli::sleep(unsigned long long milisec)
{
    QTime timer_for_sleep;

    timer_for_sleep.start();

    while (timer_for_sleep.elapsed() < milisec);
}

void cli::generate_permutation(long permutation[], bool occupied[], refer level, algorithm *alg)
{
    refer i;

    if (level == G->n)
    {
        /*for (i=0;i<G->n;i++)
        {
            printf("%u",permutation[i]);
        }*/
        alg->gcc_ccp(G, result_for_permutation, permutation);
        if (prob->count_labels(G, result_for_permutation) == 5)
        {
            for (i=0;i<G->n;i++)
            {
                sizes[i] = 0;
            }
            for (i=0;i<G->n;i++)
            {
                sizes[result_for_permutation[i]]++;
            }
            if (sizes[1] == 2 && sizes[2] == 2 && sizes[3] == 2 && sizes[4] == 2 && sizes[5] == 2 &&
                result_for_permutation[0] == result_for_permutation[5] &&
                result_for_permutation[1] == result_for_permutation[6] &&
                result_for_permutation[2] == result_for_permutation[7] &&
                result_for_permutation[3] == result_for_permutation[4] &&
                result_for_permutation[8] == result_for_permutation[9])
            {
                bad_scenario_counter++;
            }
        }
        histogram[prob->count_labels(G, result_for_permutation)]++;
        return;
    }

    for (i=0;i<G->n;i++)
    {
        if (! occupied[i])
        {
            permutation[level] = i;
            occupied[i] = true;
            generate_permutation(permutation, occupied, level + 1, alg);
            permutation[level] = 0;
            occupied[i] = false;
        }
    }
}

void cli::try_all_permutations()
{
    refer i;
    long *permutation = new long[MAX_VERTICES];
    bool *occupied    = new bool[MAX_VERTICES];
    bad_scenario_counter = 0;

    algorithm *alg = new algorithm_iggcc();
    prob = new problem_ccp();

    for (i=0;i<G->n;i++)
    {
        permutation[i] = 0;
        occupied[i] = false;
        histogram[i] = 0;
    }

    generate_permutation(permutation, occupied, 0, alg);

    delete(prob);
    delete(alg);

    delete[] permutation;
    delete[] occupied;

    printf("4: %ld, 5: %ld, 6: %ld, bad: %ld\n",histogram[4],histogram[5],histogram[6],bad_scenario_counter);
    getchar();
}

int cli::start_cli()
{
    int j;
    double max_runs = (double) (MAX_INSTANCE * MAX_ALGORITHM);
    printf("VERTEX CLIQUE COVERING PROBLEM (CCP) SOLVER\nversion 0.2\n");

    avg_brelaz = 0.0;
    avg_iggcc = 0.0;
    avg_gis = 0.0;
    avg_rls_is = 0.0;
    approx = 1.0;
    discrep = 0;

    // instance loaded from file
    choose_instance();
    compute_statistics();
    choose_algorithm();
    //try_all_permutations();


    // instance generated by an algorithm

    /*printf("w = "); scanf("%lu",&w);
    printf("n_max = "); scanf("%lu",&n_max);

    for (j=0;j<MAX_INSTANCE;j++)
    {
        generate_instance();
        compute_statistics();
        sleep(1000);
        choose_algorithm();
    }*/
    printf("Final results: IG-GCC: %0.2lf, BRE: %0.2lf, GIS: %0.2lf, RLS-IS: %0.2lf\n",avg_iggcc/max_runs,avg_brelaz/max_runs,avg_gis/max_runs,avg_rls_is/max_runs);
    printf("Approximation: %0.6lf, discrepancy: %ld\n", approx, discrep);

    return 0;
}

int cli::generate_instance()
{
    if (get_graph() != NULL)
    {
        free_graph();
    }
    generate_graph_BA_model(w,n_max);
    G = get_graph();

    return 0;
}

int cli::choose_instance()
{
    printf("File:");
    scanf("%s",&filename);

    printf("Loading the data...\n");
    if ((source = fopen(filename,"r")) == NULL)
    {
        printf("Error: An error occured during opening of the file.");
        getchar();
        return 1;
    }
    if (get_graph() != NULL)
    {
        free_graph();
    }
    input_graph(source);
    fclose(source);
    G = get_graph();

    return 0;
}

int cli::compute_statistics()
{
    refer i, maxdeg;
    unsigned long triangles;

    printf("Computing the basic statistics...\n");

    min_deg = statistics::min_degree(G);
    max_deg = statistics::max_degree(G);
    avg_deg = statistics::average_degree(G);
    stdev_deg = statistics::degree_stdev(G);
    triangles = statistics::triangles(G);

    printf("G: %s\n",filename);
    printf("|V| = %ld, |E| = %ld, d = %0.3lf\n",G->n,G->m,(double)(2*G->m)/(double)(G->n)/(double)(G->n-1));
    printf("min degree = %ld, max degree = %ld\n",min_deg,max_deg);
    printf("average degree = %0.3lf, degree stdev = %0.3lf\n",avg_deg,stdev_deg);
    printf("number of triangles = %lu\n",triangles);

    printf("degree distribution: [");
    maxdeg = statistics::degree_distribution(G, degree_distrib);
    for (i=1;i<=maxdeg;i++)
    {
        printf("%ld", degree_distrib[i]);
        if (i != maxdeg)
        {
            putchar(',');
        }
    }
    printf("]\n");

    return 0;
}

int cli::choose_algorithm()
{
    long j;
    refer current_gcc;
    refer* initial_indset = new refer[MAX_VERTICES];
    refer initial_indset_size;

    printf("Finding clique coverings using our heuristic...\n");

    timer.start();

    problem_ccp *ccp;
    refer *result;
    refer indset_size;

    // GIS
    printf("GIS:\n");
    algorithm_greedyindset *greedy_indset;
    greedy_indset = new algorithm_greedyindset();

    result = new refer[G->n];

    avg_time = 0;
    avg_iter = 0;
    succ_rate = 0;
    i = 0;
    //for (i=0;i<MAX_ALGORITHM;i++)
    {

        printf("Run %ld:",i+1);
        srand((unsigned long long) time(NULL));

        indset_size = greedy_indset->greedy_indset(G,result);

        //if (is_optimal)
        //{
        //    succ_rate++;
        //}

        printf("\nThe algorithm found a an independent set of size %u.\n",indset_size);
        avg_time += timer.elapsed() / 1000;
        avg_gis += (double) (indset_size);
        sleep(1000);
        timer.restart();


    }

    for (i=0;i<G->n;i++)
    {
        initial_indset[i] = result[i];
    }
    initial_indset_size = indset_size;

    //avg_time /= MAX_ALGORITHM;
    //avg_iter /= MAX_ALGORITHM;
    //printf("\nOn average, it took %llu seconds and %llu iterations. Success rate: %ld / 10.\n",avg_time,avg_iter,succ_rate);

    delete(greedy_indset);
    delete(result);

    // IG-GCC
    printf("IG-GCC:\n");
    algorithm_iggcc *iggcc_ccp;
    iggcc_ccp = new algorithm_iggcc();

    ccp = new problem_ccp();
    result = new refer[G->n];

    avg_time = 0;
    avg_iter = 0;
    succ_rate = 0;
    min_iggcc = MAX_VERTICES;
    for (i=0;i<MAX_ALGORITHM;i++)
    {

        printf("Run %ld:",i+1);
        srand((unsigned long long) time(NULL));

        //do
        {
            is_optimal = iggcc_ccp->iggcc_ccp(G,result,&indset_size,initial_indset,initial_indset_size);
        }
        //while (! is_optimal && timer.elapsed() < 18000000);

        if (is_optimal)
        {
            succ_rate++;
        }

        printf("\nThe algorithm needed %u cliques.\n",ccp->count_labels(G,result));
        avg_time += timer.elapsed() / 1000;
        current_gcc = ccp->count_labels(G,result);
        avg_iggcc += (double) (current_gcc);
        if (min_iggcc > current_gcc)
        {
            min_iggcc = current_gcc;
        }
        avg_rls_is += (double) (indset_size);
        if (approx < (double) current_gcc / (double) indset_size)
        {
            approx = (double) current_gcc / (double) indset_size;
            discrep = current_gcc - indset_size;
        }
        sleep(1000);
        timer.restart();

    }
    FILE *f = fopen("solution.txt", "w");
    for (i=1;i<=ccp->count_labels(G,result);i++)
    {
        for (j=0;j<G->n;j++)
        {
            if (i == result[j])
            {
                fprintf(f,"%u ",j+1);
            }
        }
        fprintf(f,"\n");
    }
    fclose(f);

    avg_time /= MAX_ALGORITHM;
    avg_iter /= MAX_ALGORITHM;
    printf("\nOn average, it took %llu seconds and %llu iterations. Success rate: %ld / %ld.\n",avg_time,avg_iter,succ_rate,MAX_ALGORITHM);

    delete(iggcc_ccp);
    delete(ccp);
    delete(result);

    delete[] initial_indset;

    // Brelazovka
    /*printf("Brelaz:\n");
    algorithm_brelaz_ccp *brelaz_algorithm;
    brelaz_algorithm = new algorithm_brelaz_ccp(G,G->n);

    ccp = new problem_ccp();
    result = new refer[G->n];

    avg_time = 0;
    avg_iter = 0;
    succ_rate = 0;
    refer mincl = G->n;
    for (i=0;i<MAX_ALGORITHM;i++)
    {

        printf("Run %ld:",i+1);
        srand((unsigned long long) time(NULL));

        brelaz_algorithm->brelaz_ccp(G,result,G->n);

        if (is_optimal)
        {
            succ_rate++;
        }

        printf("\nThe algorithm needed %u cliques.\n",ccp->count_labels(G,result));
        avg_time += timer.elapsed() / 1000;
        avg_brelaz += (double) (ccp->count_labels(G,result));
        if (mincl > ccp->count_labels(G,result))
        {;
            mincl = ccp->count_labels(G,result);
        }
        sleep(1000);
        timer.restart();

    }
    printf("minimum by Brelaz: %u\n", mincl);

    avg_time /= MAX_ALGORITHM;
    avg_iter /= MAX_ALGORITHM;
    printf("\nOn average, it took %llu seconds and %llu iterations. Success rate: %ld / 10.\n",avg_time,avg_iter,succ_rate);


    delete(brelaz_algorithm);
    delete(ccp);
    delete(result);*/

    // TABUCOL-CCP

    /*printf("k = ");
    scanf("%hu",&max_label);

    timer.start();

    // step 1 - initialiation using Brelaz's heuristic
    algorithm_brelaz_ccp *brelaz_algorithm;
    algorithm_tabucol_ccp *tabucol_ccp;
    ccp = new problem_ccp();
    result = new refer[G->n];
    brelaz_algorithm = new algorithm_brelaz_ccp(G,G->n);
    tabucol_ccp = new algorithm_tabucol_ccp();

    avg_time = 0;
    avg_iter = 0;
    succ_rate = 0;
    for (i=0;i<30;i++)
    {
        printf("Run %ld:",i+1);
        srand((unsigned long long) time(NULL));

        // Brelaz
        brelaz_algorithm->brelaz_ccp(G,result,max_label);
        ccp->label_unlabeled_vertices_randomly(G,result,max_label);

        // initialization of TabuCol-CCP
        if (i == 0)
        {
            tabucol_ccp->initialize(G,result,max_label);
        }
        else
        {
            tabucol_ccp->reinitialize(G,result,max_label);
        }

        // step 3 - optimization in TabuCol-CCP
        do
        {
            is_optimal = tabucol_ccp->tabucol_ccp(G,result,max_label,10000000/G->n);
            printf(",%lld",tabucol_ccp->get_conflicts());
        }
        while (! is_optimal && timer.elapsed() < 18000000);

        if (is_optimal)
        {
            succ_rate++;
        }

        current_time = timer.elapsed()/1000;
        avg_time += current_time;
        avg_iter += tabucol_ccp->get_iterations();
        printf("\nTook %llu seconds and %llu iterations.\n",current_time,tabucol_ccp->get_iterations());
        timer.restart();

    }

    avg_time /= 30;
    avg_iter /= 30;
    printf("\nOn average, it took %llu seconds and %llu iterations. Success rate: %ld / 30.\n",avg_time,avg_iter,succ_rate);

    tabucol_ccp->finalize();
    delete(tabucol_ccp);
    delete(brelaz_algorithm);
    delete(ccp);
    delete(result);*/

    free_graph();

    return 0;
}
