#include <stdio.h>
#include "cli.h"
#include "random_generator.h"
#include <time.h>
int main(int argc, char *argv[])
{
    /*srand((unsigned long long) time(NULL));
    random_generator gen;
    for (int i=0; i<100; i++) printf("%d ", gen.random(0,19999));
    return 0;*/
    cli *cli_instance;

    cli_instance = new cli();
    cli_instance->start_cli();
    delete(cli_instance);

    return 0;
}
