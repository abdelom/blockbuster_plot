#include "blockbuster_grid.h"
#include "blockbuster.h"
#include "sfs.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "linear_regression_f.h"


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
    int troncation;
    double *theta_flag;
    int theta_flag_count;
} Args;



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
    args->troncation = 0;
    args->theta_flag = NULL;        // <-- tableau de floats (à allouer après parsing)
    args->theta_flag_count = 0;     // <-- nombre d'éléments

    // Define long options
    static struct option long_options[] = {
        {"output_directory", required_argument, 0, 'p'},
        {"troncation", required_argument, 0, 't'},
        {"changes", required_argument, 0, 'c'},
        {"sfs", required_argument, 0, 's'},
        {"oriented", required_argument, 0, 'o'},
        {"blocks", required_argument, 0, 'b'},
        {"upper_bound", required_argument, 0, 'u'},
        {"lower_bound", required_argument, 0, 'l'},
        {"grid_size", required_argument, 0, 'n'},
        {"recent", required_argument, 0, 'r'},
        {"singleton", required_argument, 0, 'S'},
        {"theta_flag", required_argument, 0, 'T'}, // <-- nouvelle option renommée
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    // Parse command-line arguments
    int opt;
    while ((opt = getopt_long(argc, argv, "c:p:o:b:l:u:h:s:r:n:S:t:T:", long_options, NULL)) != -1)
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
        case 't':
            args->troncation = atoi(optarg);
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
            break;
        case 'T': // <-- parsing liste de float pour theta_flag
        {
            char *token = strtok(optarg, ",");
            while (token != NULL)
            {
                args->theta_flag = realloc(args->theta_flag, (args->theta_flag_count + 1) * sizeof(double));
                args->theta_flag[args->theta_flag_count] = atof(token);
                args->theta_flag_count++;
                token = strtok(NULL, ",");
            }
            break;
        }
        case 'h':
            usage(argv[0]);
            return 1;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    // Vérifications obligatoires
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

    return 0;
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


// void improve_solution(int * time_scale, solution sol, Args args, double ** sfs, int n_sample)
// {
//     int grid_size2 = args.grid_size  * GRIDREFINE;
//     double *H = generate_logarithmic_scale(grid_size2, args.upper_bound, args.lower_bound, "scenario.txt"); // Logarithmic scale
//     double **cumul_weight = cumulatve_weight_v2(n_sample, args.grid_size, H);
// }


void save_list_solution(Solution *list_solution, Solution *list_solution2, SFS sfs, char *out_file, int changes, Time_gride tg)
{
    for (int i = 0; i <= changes; i++)
    {
        if (list_solution2 == NULL)
            save_solution(list_solution[i], sfs, tg, out_file);
        else
        {
            save_solution(list_solution2[i], sfs, tg, out_file);
            clear_solution(list_solution2[i]);
        }
        clear_solution(list_solution[i]);
    }
}


int main(int argc, char *argv[])
{
    // Structure pour stocker les arguments
    Args args;
    // srand(time(0));
    // Appel de la fonction de parsing des arguments
    if (parse_args(argc, argv, &args) != 0)
        return 1; // Erreur dans le parsing
    char *outfile = construct_output_filepath(args.prefixe, "scenarios.txt");
    // Generate the time scale (either logarithmic or linear based on the `log` flag)
    SFS sfs= int_sfs(args.sfs_file, args.oriented, args.troncation, args.singleton, args.num_blocks);
    if (args.lower_bound < 0)
        args.lower_bound = 1. / (10 * sfs.n_haplotypes);
    if (args.recent > 0)
        args.lower_bound = 1. / (double)(args.recent * sfs.n_haplotypes);
    Time_gride time_grid = init_time_grid(sfs, args.grid_size, args.upper_bound, args.lower_bound);
    Flag flag = init_flag(sfs.sfs_length, args.theta_flag , args.theta_flag_count);
    clock_t start_time = clock(); // Start time measurement
    Solution *list_solution;
    if (args.theta_flag_count > 0)
    {
        list_solution = find_scenario_f(sfs, time_grid, args.changes, flag);
        save_solution(list_solution[0], sfs, time_grid, outfile);
        clear_solution(list_solution[0]);
        free(flag.thetas);
        free(flag.sfs);
        // c =  flag.n_theta - 1;
    }
    else
    {
        list_solution = find_scenario(sfs, time_grid, args.changes);
        save_list_solution(list_solution, NULL, sfs, outfile, args.changes, time_grid);
    }
    // generate_brk_combinations_f(flag.n_theta - 1, sfs, time_grid, flag);
    clock_t end_time = clock();                                                // End time measurement
    double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds
    
    printf("\n Time taken for system resolution: %f seconds\n", cpu_time_used); // Print the elapsed time
    // free_integral_grid(cumul_weight, n_sample);

    free(outfile);
    clear_time_grid(time_grid, sfs.sfs_length);
    clear_sfs(sfs);
    free(list_solution);
    return 0;
}


    // // solution *list_solution2;
    // // if (args.recent > 0)
    // //     list_solution2 = recent_infrence(list_s free(sfs);olution, args.changes, sfs, cumul_weight, size, n_sample);
    // // else
    // //     list_solution2 = NULL;
