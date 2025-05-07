#include "blockbuster_grid.h"
#include "blockbuster.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef struct
{
    int oriented;
    int singleton;
    int num_blocks;
    int changes;
    double recent;
    char *prefixe;
    char *sfs_file;
    double upper_bound;
    double lower_bound;
    int grid_size;
} Args;

/**
 * Reads the Site Frequency Spectrum (SFS) data from a file and returns it as a dynamically allocated array.
 * The function reads all floating-point numbers from the file, stores them as doubles, and returns the
 * array and its size through a pointer argument.
 *
 * @param filename The path to the file containing the SFS data, where each number represents an SFS entry.
 * @param size Pointer to an integer where the function will store the number of entries in the SFS.
 *
 * @return A dynamically allocated array of doubles containing the SFS data from the file. The caller is
 *         responsible for freeing this memory.
 *
 * Example usage:
 *     int size;
 *     double *sfs = readSFSFromFile("sfs_data.txt", &size);
 *     // Use `sfs` array here
 *     free(sfs);  // Free the allocated memory when done
 */
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

void usage(char *prog_name)
{
    fprintf(stderr, "Usage: %s -f <filename> [-o <oriented>] [-b <num_blocks>] [-u <upper_bound>] [--help]\n", prog_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -p, --prefixe_directory <output_directory>\n");
    fprintf(stderr, "  -s, --sfs <sfs>                Filename containing the Site Frequency Spectrum (SFS) (required)\n");
    fprintf(stderr, "  -o, --oriented <1|0>           Oriented (default: 0)\n");
    fprintf(stderr, "  -b, --blocks <num_blocks>      Number of blocks (default: 3)\n");
    fprintf(stderr, "  -u, --upper_bound <value>      Upper bound for time grid in Ne(0) generations (double, default: 1, must be positive)\n");
    fprintf(stderr, "  -l, --lower_bound <value>      lower bound for time grid in Ne(0) generations (double, default: 1e-4, must be positive)\n");
    fprintf(stderr, "  -n, --grid_size <value>        number of tme point between lower bound and uppper bount(integer, default: 33\n");
    fprintf(stderr, "  -c, --changes <value>      \n");
    fprintf(stderr, "      --help                     Display this help and exit\n");
}

int parse_args(int argc, char *argv[], Args *args)
{
    // Default values
    args->oriented = 0;
    args->num_blocks = 1;
    args->prefixe = NULL;
    args->sfs_file = NULL;
    args->upper_bound = 1.0;
    args->lower_bound = -1.0;
    args->changes = 5;
    args->recent = -1.;
    args->grid_size = 35;
    args->singleton = 1;

    // Define long options
    static struct option long_options[] = {
        {"output_directory", required_argument, 0, 'p'},
        {"changes", required_argument, 0, 'c'},
        {"sfs", required_argument, 0, 's'},
        {"oriented", required_argument, 0, 'o'},
        {"blocks", required_argument, 0, 'b'},
        {"upper_bound", required_argument, 0, 'u'},
        {"lower_bound", required_argument, 0, 'l'},
        {"grid_size", required_argument, 0, 'n'},
        {"recent", required_argument, 0, 'r'},
        {"help", no_argument, 0, 'h'},
        {"singleton", required_argument, 0, 'S'},
        {0, 0, 0, 0} // End of list
    };

    // Parse command-line arguments
    int opt;
    while ((opt = getopt_long(argc, argv, "c:p:o:b:l:u:h:s:r:n:S:", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'p':
            args->prefixe = optarg;
            break;
        case 'n':
            args->grid_size = atoi(optarg);
            break;
        case 's':
            args->sfs_file = optarg;
            break;
        case 'o':
            args->oriented = atoi(optarg);
            break;
        case 'S':
            args->singleton = atoi(optarg);
            break;
        case 'r':
            args->recent = atof(optarg);
            break;
        case 'b':
            args->num_blocks = atoi(optarg);
            break;
        case 'c':
            args->changes = atoi(optarg);
            if (args->changes <= 0)
            {
                fprintf(stderr, "Error:\n");
                usage(argv[0]);
                return 1;
            }
            break;
        case 'u':
            args->upper_bound = atof(optarg);
            if (args->upper_bound <= 0)
            {
                fprintf(stderr, "Error: upper bound (-u or --upper-bound) must be positive.\n");
                usage(argv[0]);
                return 1;
            }
            break;
        case 'l':
            args->lower_bound = atof(optarg);
            // if (args->lower_bound <= 0)
            // {
            //     fprintf(stderr, "Error: lower bound (-l or --lower-bound) must be positive.\n");
            //     usage(argv[0]);
            //     return 1;
            // }
            break;
        case 'h':
            usage(argv[0]);
            return 1;
        default: /* '?' */
            usage(argv[0]);
            return 1;
        }
    }

    // Check that filename is provided
    if (args->prefixe == NULL)
    {
        fprintf(stderr, "Error: prefixe is required.\n");
        usage(argv[0]);
        return 1;
    }

    if (args->sfs_file == NULL)
    {
        fprintf(stderr, "Error: sfs is required.\n");
        usage(argv[0]);
        return 1;
    }

    return 0; // No errors
}

int ensure_directory_exists(const char *path)
{
    struct stat st = {0};

    // Check if directory exists
    if (stat(path, &st) == -1)
    {
        // Directory does not exist, so create it
        if (mkdir(path, 0700) == -1)
        { // 0700 gives rwx permissions to the owner
            perror("Error creating directory");
            return -1; // Return an error code if mkdir fails
        }
    }
    return 0; // Success
}

char *construct_output_filepath(const char *prefix, const char *filename)
{
    if (ensure_directory_exists(prefix) == -1)
    {
        fprintf(stderr, "Failed to ensure directory exists: %s\n", prefix);
        return NULL;
    }

    // Calculer la taille nécessaire pour le chemin complet
    size_t len_prefix = strlen(prefix);
    size_t len_filename = strlen(filename);
    size_t total_len = len_prefix + len_filename + 2; // +1 pour '/' et +1 pour le terminateur NULL

    // Allouer de la mémoire pour la chaîne complète
    char *out_file = malloc(total_len * sizeof(char));
    if (out_file == NULL)
    {
        perror("Error allocating memory");
        return NULL;
    }

    // Construire le chemin complet
    snprintf(out_file, total_len, "%s/%s", prefix, filename);

    return out_file;
}



void improve_solution(int * time_scale, solution sol, Args args, double ** sfs, int n_sample)
{
    int grid_size2 = args.grid_size  * GRIDREFINE;
    double *H = generate_logarithmic_scale(grid_size2, args.upper_bound, args.lower_bound, "scenario.txt"); // Logarithmic scale
    double **cumul_weight = cumulatve_weight_v2(n_sample, args.grid_size, H);
}


void save_list_solution(solution *list_solution, solution *list_solution2, int n_sample, int sfs_length, char *out_file, double **sfs, int changes)
{
    double const_ren = (sfs[0][0] + sfs[1][0]) / sfs[0][0];
    if ((int)const_ren == 2)
        const_ren = 1.;
    for (int i = 0; i <= changes; i++)
    {
        if (list_solution2 == NULL)
            save_solution(list_solution[i], n_sample, out_file, const_ren, sfs_length);
        else
        {
            save_solution(list_solution2[i], n_sample, out_file, const_ren, sfs_length);
            // clear_solution(list_solution2[i]);
        }
        clear_solution(list_solution[i]);
    }
}

int main(int argc, char *argv[])
{
    // Structure pour stocker les arguments
    Args args;

    // Appel de la fonction de parsing des arguments
    if (parse_args(argc, argv, &args) != 0)
        return 1; // Erreur dans le parsing

    int size;
    // Lecture du fichier et application des opérations sur SFS
    double **sfs = malloc(sizeof(double *) * 2);
    sfs[0] = readSFSFromFile(args.sfs_file, &size);
    int n_sample = size + 1; // n_sample doit être spécifié
    if (args.lower_bound < 0)
        args.lower_bound = 1. / (10 * n_sample);
    if (args.recent > 0)
        args.lower_bound = 1. / (double)(args.recent * n_sample);
    char *outfile = construct_output_filepath(args.prefixe, "scenarios.txt");
    // Generate the time scale (either logarithmic or linear based on the `log` flag)
    double *H = generate_logarithmic_scale(args.grid_size * GRIDREFINE , args.upper_bound, args.lower_bound, outfile); // Logarithmic scale
    double **cumul_weight = cumulatve_weight_v2(n_sample, args.grid_size * GRIDREFINE, H);
    // if(!args.singleton)
    //     sigleton_ignore(sfs, cumul_weight, args.grid_size);
    if(!args.oriented)
    {
        fold_sfs(sfs, cumul_weight, size, args.grid_size);
        size = size / 2 + size % 2;
    }
    // else{
    //     if(!args.singleton){
    //         singleton_erased(sfs, cumul_weight, size);
    //         size -= 1;
    //     }
    // }
    sfs[1] = test_split(sfs[0], size, args.num_blocks);
    for (int i = 0; i < size; i++)
        printf("%f %f \n", sfs[0][i], sfs[1][i]);
    clock_t start_time = clock(); // Start time measurement

    solution *list_solution = find_scenario(size, cumul_weight, sfs, args.grid_size, n_sample, args.changes);
    solution *list_solution2;
    if (args.recent > 0)
        list_solution2 = recent_infrence(list_solution, args.changes, sfs, cumul_weight, size, n_sample);
    else
        list_solution2 = NULL;
    clock_t end_time = clock();                                                // End time measurement
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds

    save_list_solution(list_solution, list_solution2, n_sample, size, outfile, sfs, args.changes);
    printf("Time taken for system resolution: %f seconds\n", cpu_time_used); // Print the elapsed time
    free_integral_grid(cumul_weight, n_sample);
    free(sfs[1]);
    free(sfs[0]);
    free(sfs);
    free(outfile);
    free(list_solution);
    return 0;
}
