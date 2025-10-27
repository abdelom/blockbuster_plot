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
    int epochs;             
    double recent;          
    char *prefixe;          
    char *sfs_file;         
    double upper_bound;     
    double lower_bound;     
    int grid_size;          
    int troncation;         
    double *parameters_flag; 
    int parameters_count; 
    double genome_length;
    double mutation_rate;
    double generation_time;   
} Args;


int parse_parameters_flag(const char *str, double **out_array, int *out_count) {
    // Vérifier chaîne vide
    if (str == NULL || *str == '\0') {
        return 0;  // Rien à faire, on retourne 0 sans modifier les autres arguments
    }

    char *copy = strdup(str);
    if (!copy) return 1;

    char *pipe_pos = strchr(copy, '|');
    if (!pipe_pos) {
        fprintf(stderr, "Error: parameters_flag must contain a '|'\n");
        free(copy);
        return 1;
    }

    *pipe_pos = '\0';
    char *left = copy;
    char *right = pipe_pos + 1;

    int left_count = 0, right_count = 0;
    char *token;

    char *left_copy = strdup(left);
    char *right_copy = strdup(right);
    if (!left_copy || !right_copy) {
        perror("strdup");
        free(copy);
        free(left_copy);
        free(right_copy);
        return 1;
    }

    // Compter éléments gauche
    token = strtok(left_copy, ",");
    while (token != NULL) {
        left_count++;
        token = strtok(NULL, ",");
    }

    // Compter éléments droite
    token = strtok(right_copy, ",");
    while (token != NULL) {
        right_count++;
        token = strtok(NULL, ",");
    }

    free(left_copy);
    free(right_copy);

    if (left_count != right_count + 1) {
        fprintf(stderr, "Error: Left side must have right_count + 1 elements (%d + 1 = %d), but got %d\n",
                right_count, right_count + 1, left_count);
        free(copy);
        return 1;
    }

    // Allouer tableau (concaténer gauche + droite)
    *out_count = left_count + right_count;
    *out_array = malloc((*out_count) * sizeof(double));
    if (!*out_array) {
        perror("malloc");
        free(copy);
        return 1;
    }

    // Remplir tableau
    int idx = 0;
    token = strtok(left, ",");
    while (token != NULL) {
        (*out_array)[idx++] = atof(token);
        token = strtok(NULL, ",");
    }
    token = strtok(right, ",");
    while (token != NULL) {
        (*out_array)[idx++] = atof(token);
        token = strtok(NULL, ",");
    }

    free(copy);
    return 0;
}


// Affichage de l'aide
void usage(char *prog_name)
{
    fprintf(stderr, "Usage: %s -f <filename> [-o <oriented>] [-b <num_blocks>] [-u <upper_bound>] [--help]\n", prog_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -p, --prefixe_directory <output_directory>\n");
    fprintf(stderr, "  -s, --sfs <sfs>                Filename containing the Site Frequency Spectrum (SFS) (required)\n");
    fprintf(stderr, "  -o, --oriented <1|0>           Oriented (default: 0)\n");
    fprintf(stderr, "  -b, --blocks <num_blocks>      Number of blocks (default: 3)\n");
    fprintf(stderr, "  -u, --upper_bound <value>      Upper bound for time grid in Ne(0) generations (default: 1, must be positive)\n");
    fprintf(stderr, "  -l, --lower_bound <value>      Lower bound for time grid (default: 1e-4)\n");
    fprintf(stderr, "  -n, --grid_size <value>        Number of time points between lower and upper bound (default: 35)\n");
    fprintf(stderr, "  -e, --epochs_max <value>       Maximum number of epochs (default: 5)\n");
    fprintf(stderr, "  -r, --recent <value>           Recent time (optional)\n");
    fprintf(stderr, "  -S, --singleton <1|0>          Include singletons (default: 1)\n");
    fprintf(stderr, "  -t, --troncation <0|1>         Troncation option\n");
    fprintf(stderr, "  -P, --parameters_flag <string>  Parameters formatted as '0,0|0,0,0'\n");
    fprintf(stderr, "  -m, --mutation_rate         Mutation rate per generation per base pair (Time scale in unit of Ne genration if not specifie)\n");
    fprintf(stderr, "  -L, --Genome_length         Genome length in base pair (Time scale in unit of Ne genration if not specifie)\n");
    fprintf(stderr, "  -g, --Generation_time         generation time in year (Time scale in unit of generation if not specifie)\n");
    fprintf(stderr, "      --help                     Display this help and exit\n");
}

int parse_args(int argc, char *argv[], Args *args)
{
    // Valeurs par défaut
    args->oriented = 0;
    args->num_blocks = 1;
    args->prefixe = NULL;
    args->sfs_file = NULL;
    args->upper_bound = 1.0;
    args->lower_bound = 1e-4;
    args->epochs = 5;
    args->recent = -1.0;
    args->grid_size = 35;
    args->singleton = 1;
    args->troncation = 0;
    args->parameters_flag = NULL;
    args->parameters_count = 0;
    args->genome_length = -1.; // Default genome lengt  
    args->mutation_rate = -1.; // Default mutation rate
    args->generation_time = -1.; // Default generation time in years

    // Options longues
    static struct option long_options[] = {
        {"output_directory", required_argument, 0, 'p'},
        {"troncation", required_argument, 0, 't'},
        {"epochs", required_argument, 0, 'e'},
        {"sfs", required_argument, 0, 's'},
        {"oriented", required_argument, 0, 'o'},
        {"blocks", required_argument, 0, 'b'},
        {"upper_bound", required_argument, 0, 'u'},
        {"lower_bound", required_argument, 0, 'l'},
        {"grid_size", required_argument, 0, 'n'},
        {"mutation_rate", required_argument, 0, 'm'},
        {"genome_length", required_argument, 0, 'L'},
        {"generation_time", required_argument, 0, 'g'},
        {"recent", required_argument, 0, 'r'},
        {"singleton", required_argument, 0, 'S'},
        {"parameters_flag", required_argument, 0, 'P'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "e:p:o:b:l:u:h:s:r:n:S:t:P:m:L:g:", long_options, NULL)) != -1)
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

            case 'm':
                args->mutation_rate = atof(optarg);
                break;
            case 'L':
                args->genome_length = atof(optarg);
                break;
            case 'g':
                args->generation_time = atof(optarg);
                break;
            case 'e':
                args->epochs = atoi(optarg);
                if (args->epochs <= 0) {
                    fprintf(stderr, "Error:\n");
                    usage(argv[0]);
                    return 1;
                }
                break;
            case 'u':
                args->upper_bound = atof(optarg);
                if (args->upper_bound <= 0) {
                    fprintf(stderr, "Error: upper bound must be positive.\n");
                    usage(argv[0]);
                    return 1;
                }
                break;
            case 'l':
                args->lower_bound = atof(optarg);
                break;
            case 'P':
                if (parse_parameters_flag(optarg, &args->parameters_flag, &args->parameters_count) != 0) {
                    return 1;
                }
                break;
            case 'h':
                usage(argv[0]);
                return 1;
            default:
                usage(argv[0]);
                return 1;
        }
    }

    // Vérifications obligatoires
    if (args->prefixe == NULL) {
        fprintf(stderr, "Error: prefixe is required.\n");
        usage(argv[0]);
        return 1;
    }

    if (args->sfs_file == NULL) {
        fprintf(stderr, "Error: sfs file is required.\n");
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


void save_list_solution(Solution *list_solution, SFS sfs, char *out_file, int epochs, Time_gride tg, Args args)
{
    for (int i = 0; i < epochs; i++)
    {
        // if (list_solution2 == NULL)
        save_solution(list_solution[i], sfs, tg, out_file, args.mutation_rate, args.generation_time, args.genome_length);
        // else
        // {
        //     save_solution(list_solution2[i], sfs, tg, out_file);
        //     clear_solution(list_solution2[i]);
        // }
        clear_solution(list_solution[i]);
    }
}


double * times_c(Solution s)
{
    double *H = malloc(s.nb_breakpoints * sizeof(double));
    double time = 0.0, lb = 0.0, delta_time = 0.0;
    for (int i = 0; i < s.nb_breakpoints; i ++)
    {
        delta_time = (s.time[i] - lb) * s.thetas[0] / s.thetas[i];
        time += delta_time;
        // printf("\n%faa\n", s.time[i]);
        lb = s.time[i];
        H[i] = time;
    }
    return H;
}

static double normal01(void) {
    static int has_spare = 0;
    static double spare;
    if (has_spare) {
        has_spare = 0;
        return spare;
    } else {
        double u, v, s;
        do {
            u = 2.0 * drand48() - 1.0;
            v = 2.0 * drand48() - 1.0;
            s = u*u + v*v;
        } while (s == 0.0 || s >= 1.0);
        double mul = sqrt(-2.0 * log(s) / s);
        spare = v * mul;
        has_spare = 1;
        return u * mul;
    }
}

static double clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

Solution metropolis_hastings(Solution sol_init, SFS sfs, Time_gride tg, int iterations)
{
    const double log_proposal_sd = 0.01;
    const double H_min = 1e-4;
    const double H_max = 2.0;
    double log_max = -INFINITY;
    static int seeded = 0;
    if (!seeded) {
        srand48((unsigned) time(NULL) ^ (unsigned) getpid());
        seeded = 1;
    }

    Solution current = copy_solution(sol_init);
    double *H = times_c(sol_init);
    int len = current.nb_breakpoints;

    Solution *chaine = malloc(iterations * sizeof(Solution));
    double *H_new = malloc(len * sizeof(double));
    Time_gride tg2 = init_time_grid_H_wik(tg, sfs, H);
    memcpy(H_new, H, len * sizeof(double));

    for (int i = 0; i < iterations; ++i) {

        // ---- Proposer de nouvelles valeurs pour TOUS les H ----
        double *H_prop = malloc(len * sizeof(double));
        for (int k = 0; k < len; ++k) {
            double logH = log(H_new[k]);
            double prop = logH + normal01() * log_proposal_sd;
            double Hcand = exp(prop);
            Hcand = clamp(Hcand, H_min, H_max);
            H_prop[k] = Hcand;
        }

        // ---- Calcul du candidat ----
        tg2.time_scale = H_prop;
        cumulatve_weight_v2(sfs.n_haplotypes, sfs.n_haplotypes - 1, current.nb_breakpoints, &tg2);

        Solution cand = init_solution_size(current.nb_breakpoints, 1);
        cand.breakpoints[current.nb_breakpoints] = current.nb_breakpoints + 1;
        system_resolution(&cand, sfs, tg2);

        double loglike_old = current.log_likelihood;
        double loglike_new = cand.log_likelihood;
        double log_accept_ratio = loglike_new - loglike_old;

        // ---- Acceptation / Rejet ----
        int accept = 0;
        if (log_accept_ratio >= 0.0 || log(drand48()) < log_accept_ratio)
            accept = 1;

        if (accept) {
            current = copy_solution(cand);
            memcpy(H_new, H_prop, len * sizeof(double));
        } 

       

        if (i %1000 == 0){
            chaine[i] = copy_solution(current);
            printf(" \n%.6f %f %f %f, ", current.thetas[0],  current.thetas[1], current.thetas[2], current.log_likelihood);
            for(int l =0; l < current.nb_breakpoints; l++){
                printf("%.6f ", H_prop[l]);
            }
        }
         free(H_prop);
    }

    printf("Final logL = %.6f\n", log_max);
    return current;
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
    clock_t start_time = clock(); // Start time measurement
    Solution *list_solution;
    if (args.parameters_count > 0)
    {
        Flag flag = init_flag(sfs.sfs_length, args.parameters_flag , args.parameters_count);
        list_solution = find_scenario_f(sfs, time_grid, flag);
        save_solution(list_solution[0], sfs, time_grid, outfile, args.mutation_rate, args.generation_time, args.genome_length);
        clear_solution(list_solution[0]);
        free(flag.thetas);
        free(flag.sfs);
        // c =  flag.n_theta - 1;
    }
    else
    {
        list_solution = find_scenario(sfs, time_grid, args.epochs);
        // clock_t start_time = clock(); // Start time measurement
        // metropolis_hastings(list_solution[6], sfs, time_grid, 5e5);
        // clock_t end_time = clock();                                                // End time measurement
        // double cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC; // Calculate the time in seconds
        // printf("\ns Time taken for system resolution: %f seconds\n", cpu_time_used); // Print the elapsed time
        save_list_solution(list_solution, sfs, outfile, args.epochs, time_grid, args);
    }
    // generate_brk_combinations_f(flag.n_theta - 1, sfs, time_grid, flag);
    // free_integral_grid(cumul_weight, n_sample);

    free(outfile);
    clear_time_grid(time_grid, sfs.sfs_length);
    clear_sfs(sfs);
    free(list_solution);
    free(time_grid.wik);
    return 0;
}


    // // solution *list_solution2;
    // // if (args.recent > 0)
    // //     list_solution2 = recent_infrence(list_s free(sfs);olution, args.epochs, sfs, cumul_weight, size, n_sample);
    // // else
    // //     list_solution2 = NULL;
