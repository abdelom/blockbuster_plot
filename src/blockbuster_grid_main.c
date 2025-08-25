
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "blockbuster_grid.h"
#include <getopt.h>

/**
 * Prints the usage information for the program, detailing all the available command-line options.
 */
void print_usage()
{
    printf("Usage: program -n n_sample [-u upper_bound] [-l lower_bound] [-g grid_size] [-s scale] [-o outputfile]\n");
    printf("  -n n_sample       : Required. The sample size (positive integer).\n");
    printf("  -u upper_bound    : Optional. The upper bound for the time scale (default = 1.0).\n");
    printf("  -l lower_bound    : Optional. The lower bound for the time scale (default = 1e-4).\n");
    printf("  -g grid_size      : Optional. The size of the grid representing time intervals (default = 35).\n");
    printf("  -s scale          : Optional. Time scale type: 1 for logarithmic, 0 for linear (default = 1).\n");
    printf("  -o outputfile     : Optional. The output file name for saving results (default = 'grid.txt').\n");
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
    double *H;
    if(scale)
        H = generate_logarithmic_scale(grid_size, upper_bound, lower_bound, outputfile);
    else
        H = generate_linear_scale(grid_size, upper_bound, lower_bound, outputfile);
    // double ** weight_grid = cumulatve_weight_v2(n_sample, grid_size, H);
    Time_gride tg = init_time_grid_H(n_sample, grid_size, H);
    save_cumulated_weight(n_sample, grid_size + 2, tg.cumulative_bl, outputfile);
    clear_time_grid(tg, n_sample - 1);
    return 0;
}