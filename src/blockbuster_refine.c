#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "blockbuster.h"
#include "blockbuster_grid.h"
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

typedef struct
{
    int oriented;
    char *filin;
    char *sfs_file;
} Args;


// Function to initialize a solution with a specified number of breakpoints
solution init_solution_size_brk(int nb_breakpoints, int*brk, int grid_size)
{
    printf("%d aa\n", nb_breakpoints);
    solution sol;
    sol.nb_breakpoints = nb_breakpoints;
    // Allocate memory for `nb_breakpoints + 1` to include initial and final states or boundaries, the last one corespond to the infinity
    sol.breakpoints = realloc(brk, (nb_breakpoints + 2) * sizeof(int));
    sol.thetas = calloc(sizeof(double), (nb_breakpoints + 1)); // Allocate space for each interval’s mutation rate
    sol.se_thetas = NULL;
    sol.residues = NULL;
    sol.breakpoints = brk; // Initialize breakpoints as a simple sequence; this may be customized to represent actual change times
    sol.breakpoints[nb_breakpoints] = grid_size + 1;
    return sol;
}


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
    free_integral_grid(cumul_weight, n_sample);
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


int commence_par_breakpoint(const char *ligne) {
    return strncmp(ligne, " breakpoint", strlen(" breakpoint")) == 0;
}

int extraire_breakpoints(const char *ligne, int **breakpoints) {
    const char *mot_cle = " breakpoints:";
    char *position = strstr(ligne, mot_cle);
    if (!position) {
        *breakpoints = NULL;
        return 0;
    }

    position += strlen(mot_cle); // Aller après "breakpoints:"
    int *resultat = malloc(sizeof(int) * 100); // capacité initiale
    int count = 0;

    while (*position) {
        while (*position && !isdigit(*position) && *position != '-') {
            position++;
        }
        if (*position) {
            int valeur;
            if (sscanf(position, "%d", &valeur) == 1) {
                resultat[count++] = valeur;
            }
            while (*position && (isdigit(*position) || *position == '-')) {
                position++;
            }
        }
    }

    *breakpoints = realloc(resultat, count * sizeof(int)); // ajuster la taille
    return count;
}

void parse_scenario(char * scenarios)
{
    FILE *file = fopen(scenarios, "r");

    char line[100000]; // Taille maximale d'une ligne, ajustez selon vos besoins
    double * time = malloc(sizeof(double));
    int grid_size = 0, flag = 1;
    solution sol;
    while (fgets(line, sizeof(line), file)) {
        // Traiter la ligne
       
        if(atof(line) == 0.) flag = 0;
        if(flag)
        {
        time[grid_size] = (double)atof(line);
        grid_size ++;
        time = realloc(time, (grid_size + 1) * sizeof(double));
    // printf("%lf \n", time[i - 1]);
        }
        else
        {
            if(commence_par_breakpoint(line))
            {
                int *breaks = NULL;
                int nb_breakpoint = extraire_breakpoints(line, &breaks);
                if(nb_breakpoint != 0)
                {solution sol  = init_solution_size_brk(nb_breakpoint, breaks, grid_size);
                clear_solution(sol);}
            }
        }
         
    }
    free(time);
    // Fermer le fichier
    fclose(file);
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
     clear_solution(sol);
    parse_scenario("test/scenarios.txt");
    // Lire le fichier ligne par ligne


    const char *ligne = " breakpoints:  ";
    int *breaks = NULL;
    int nb = extraire_breakpoints(ligne, &breaks);

    printf("Nombre d'entiers : %d\n", nb);
    for (int i = 0; i < nb; ++i) {
        printf("%d ", breaks[i]);
    }
    printf("\n");

    free(breaks);
    free(sfs[0]);
    free(sfs[1]);
    free(sfs);

    return 0;
    return 0;
}
