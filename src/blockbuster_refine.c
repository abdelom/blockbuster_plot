#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "blockbuster.h"
#include "blockbuster_grid.h"

#include <ctype.h>

typedef struct
{
    int oriented;
    char *filin;
    char *sfs_file;
} Args;

int parse_args(int argc, char *argv[], Args *args)
{
    // Default values
    args->oriented = 0;
    args->sfs_file = NULL;
    args->filin = NULL;

    // Define long options
    static struct option long_options[] = {
        {"sfs", required_argument, 0, 's'},
        {"oriented", required_argument, 0, 'o'},
        {"filin", required_argument, 0, 'f'},
        {0, 0, 0, 0} // End of list
    };

    // Parse command-line arguments
    int opt;
    while ((opt = getopt_long(argc, argv, "f:o:h:s:", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'f':
            args->filin = optarg;
            break;
        case 's':
            args->sfs_file = optarg;
            break;
        case 'o':
            args->oriented = atoi(optarg);
            break;



        case 'h':
            printf("help is coming \n");
            return 1;
        default: /* '?' */
            printf("bad usage \n");
            return 1;
        }
    }
    // Check that filename is provided
    if (args->filin == NULL)
    {
        fprintf(stderr, "file is required is required.\n");
        return 1;
    }

    if (args->sfs_file == NULL)
    {
        fprintf(stderr, "Error: sfs is required.\n");
        return 1;
    }

    return 0; // No errors
}


solution solution_from_times(double * times, double ** sfs, int n_sample, int oriented, int grid_size)
{
    clock_t start_time = clock(); // Start time measurement
    double **cumul_weight = cumulatve_weight_v2(n_sample, grid_size, times);
    clock_t end_time = clock();                                                // End time measurement
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds
    printf("Time taken for grid: %f seconds\n", cpu_time_used); // Print the elapsed time
     start_time = clock(); // Start time measurement
    int sfs_length = n_sample - 1;
    solution sol = init_solution_size(grid_size);
    sol.breakpoints[grid_size] = grid_size + 1;
    if(!oriented)
    {
        fold_sfs(sfs, cumul_weight, sfs_length, grid_size);
        sfs_length = sfs_length / 2 + sfs_length % 2;
    }
   
    system_resolution(&sol, sfs, cumul_weight, sfs_length);
    end_time = clock();                                                // End time measurement
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds
    printf("Time taken for system resolution: %f seconds\n", cpu_time_used); // Print the elapsed time
    return sol;
}

double *readSFSFromFile(const char *filename, int *size)
{
    // Open the file in read mode
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1); // Exit if file cannot be opened
    }
    // Count the number of floating-point numbers in the file
    float number;
    int count = 0;
    while (fscanf(file, "%f", &number) == 1)
        count++;
    // Reset file pointer to the beginning of the file to read the values into the array
    fseek(file, 0, SEEK_SET);
    // Allocate memory for the SFS array based on the counted number of entries
    double *sfs = (double *)malloc(count * sizeof(double));
    if (sfs == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1); // Exit if memory allocation fails
    }
    // Read each float from the file into the array as a double
    int i = 0;
    while (fscanf(file, "%lf", &sfs[i]) == 1)
        i++;
    *size = count; // Set the size to the number of entries in the array
    // Close the file and return the allocated array
    fclose(file);
    return sfs;
}


// Programme principal
int main(int argc, char * argv[]) {
    
    Args args;
    double * time = malloc(sizeof(double));
    time[0] = 0.1;
    if (parse_args(argc, argv, &args) != 0)
        return 1; // Erreur dans le parsing
    int size;
    // Lecture du fichier et application des opérations sur SFS
    double **sfs = malloc(sizeof(double *) * 2);
    sfs[0] = readSFSFromFile(args.sfs_file, &size);
    sfs[1] = test_split(sfs[0], size, 1);
    int n_sample = size + 1; // n_sample doit être spécifié
    
    clock_t start_time = clock(); // Start time measurement
    solution sol = solution_from_times(time, sfs,n_sample, 1, 1);
    clock_t end_time = clock();                                                // End time measurement
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds
    printf("Time taken for system resolution: %f seconds\n", cpu_time_used); // Print the elapsed time
    printf("%s\n", args.filin);

    return 0;
}
