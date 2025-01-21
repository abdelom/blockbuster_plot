
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "blockbuster_grid.h"
#include <getopt.h>

void print_usage()
{
    printf("Usage: program -n n_sample [-u upper_bound] [-g grid_size] [-c num_core]\n");
    printf("  -n n_sample       : Required. The sample size (positive integer).\n");
    printf("  -u upper_bound    : Optional. The upper bound (default = 1.0).\n");
}

int main(int argc, char *argv[])
{
    // Default values
    int n_sample = -1; // n_sample must be specified
    double upper_bound = 1.0;
    double lower_bound = 1e-4;
    int grid_size = 35;
    int num_threads = 1; // Set default to the number of available CPU cores!
    int scale = 1;
    char *outputfile = "grid.txt";
    int opt;
    while ((opt = getopt(argc, argv, "n:u:o:l:g:s:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            n_sample = atoi(optarg);
            break;
        case 'u':
            upper_bound = atof(optarg);
            break;
        case 'o':
            outputfile = optarg; // Set the number of threads
            break;
         case 'l':
            lower_bound = atof(optarg); // Set the number of threads
            break;
        case 'g':
            grid_size = atoi(optarg); // Set the number of threads
            break;
        case 's':
            scale = atoi(optarg); // Set the number of threads
            break;
        default:
            print_usage();
            exit(EXIT_FAILURE);
        }
    }
    if (n_sample <= 0)
    {
        printf("Error: The sample size (n_sample) is required and must be a positive integer.\n");
        print_usage();
        exit(EXIT_FAILURE);
    }
    double ** weight_grid = cumulatve_weight_v2(n_sample, grid_size, upper_bound, lower_bound, "o", scale);
    save_cumulated_weight(n_sample, grid_size + 2, weight_grid, outputfile);
    free_integral_grid(weight_grid, n_sample);
    printf("%d\n", grid_size);

    return 0;
}