#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "blockbuster_grid.h"

void print_usage(char *program_name, FILE *std)
{
    fprintf(std, "Usage: %s --n_samples <n_samples> --l_genome <l_genome> --theta <theta> --oriented <oriented> --n_sim <n_sim> --model <model> --random_seed <random_seed> --parameters <parameter1,parameter2> --outputfile <outputfile>\n\n", program_name);
    fprintf(std, "Options:\n");
    fprintf(std, "  -n, --n_samples <n_samples>: Number of sampled haploides individuals\n");
    fprintf(std, "  -l, --l_genome <l_genome>: Length of the genome\n");
    fprintf(std, "  -t, --theta <theta>: Theta/population mutatiolnal rate per sites value\n");
    fprintf(std, "  -f, --oriented <oriented>: Optional, Orientation flag (0 for unoriented, 1 for oriented), if not specified set by default to 1\n");
    fprintf(std, "  -r, --random_seed <random_seed>: Optional, Random seed to init the pseudo random generator, if not specified the random seed is set to the current time\n");
    fprintf(std, "  -p, --parameters <parameter1,parameter2>: Model parameters, if constant_piecewise is specified then the number of paramter must be even times of changes comming before intensities of changes. if linear or exponential are specified then a single paramter is expecteed\n");
    fprintf(std, "  -o, --outputfile <outputfile>: Output file name in wich the sfs will be right\n");
}

void trim(char *str)
{
    // Trim leading spaces and newlines
    char *firstNonSpace = str;
    while (*firstNonSpace && isspace(*firstNonSpace))
        firstNonSpace++;
    memmove(str, firstNonSpace, strlen(firstNonSpace) + 1);
    // Trim trailing spaces and newlines
    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end))
        end--;
    // Null-terminate the trimmed string
    *(end + 1) = '\0';
}

typedef struct
{
    int n_samples;
    float l_genome;
    double theta;
    int oriented;
    int seed;
    char outfile[100];
    double *parameters;
    int n_parameters;
    int noised;
    int header;
} Args;

void init_sfs_from_args(int argc, char *argv[], Args *args)
{
    int opt;
    args->n_samples = -1;
    args->theta = -1.0;
    args->oriented = 1;
    args->seed = -1;
    args->n_parameters = 0;
    args->noised = 1;
    args->parameters = NULL;
    args->header = 0;
    args->outfile[0] = '\0'; // Initialize outfile to an empty string

    static struct option long_options[] = {
        {"n_samples", required_argument, 0, 'n'},
        {"theta", required_argument, 0, 't'},
        {"oriented", required_argument, 0, 'f'},
        {"random_seed", required_argument, 0, 'r'},
        {"outputfile", required_argument, 0, 'o'},
        {"parameters", required_argument, 0, 'p'},
        {"noised", required_argument, 0, 'e'},
        {"header", no_argument, 0, 'H'},
        {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "n:t:f:o:p:r:e:H", long_options, NULL)) != -1)
    {
        switch (opt)
        {
        case 'n':
            args->n_samples = atoi(optarg);
            break;
        case 't':
            args->theta = atof(optarg);
            break;
        case 'f':
            args->oriented = atoi(optarg);
            break;
        case 'r':
            args->seed = optarg ? atoi(optarg) : -1;
            break;
        case 'o':
            snprintf(args->outfile, sizeof(args->outfile), "%s", optarg);
            break;
        case 'e':
            args->noised = atoi(optarg);
            break;
        case 'p':
            args->parameters = malloc(100 * sizeof(double)); // Initial allocation
            if (!args->parameters)
            {
                perror("Memory allocation failed for parameters");
                exit(EXIT_FAILURE);
            }
            char *param = strtok(optarg, ",");
            while (param != NULL)
            {
                args->parameters[args->n_parameters++] = atof(param);
                param = strtok(NULL, ",");
            }
            break;
        case 'H':
            args->header = 1;
            break;
        default:
            print_usage(argv[0], stderr);
            exit(EXIT_FAILURE);
        }
    }

    if (args->n_samples == -1 || args->n_parameters == 0)
    {
        printf("Incomplete configuration\n");
        print_usage(argv[0], stderr);
        exit(EXIT_FAILURE);
    }
}

void free_args(Args *args)
{
    if (args->parameters)
    {
        free(args->parameters);
    }
}

void write_sfs_to_file(double *sfs, int n_samples, int oriented, const char *outfile)
{
    FILE *file = fopen(outfile, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: Unable to open output file.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n_samples - 1; i++)
    {
        if (oriented)
        {
            fprintf(file, "%.9f ", sfs[i]);
        }
        else
        {
            if (i < n_samples / 2)
            {
                fprintf(file, "%.9f ", sfs[i]);
            }
            else
            {
                fprintf(file, "%d ", 0);
            }
        }
    }

    fclose(file);
}

double *absolute_to_relative_times(double *parameters, int n_parameters)
{
    double *H = malloc(sizeof(double) * n_parameters / 2);
    H[0] = parameters[0];

    for (int i = 1; i < n_parameters / 2; i++)
    {
        H[i] = H[i - 1] + (parameters[i] - parameters[i - 1]) / parameters[n_parameters / 2 + i - 1];
    }
    return H;
}

double generate_normal_random(double mu, double sigma)
{
    double u1 = rand() / (double)RAND_MAX; // Génère un nombre aléatoire entre 0 et 1
    double u2 = rand() / (double)RAND_MAX;
    // Génère un autre nombre aléatoire entre 0 et 1
    double z = sqrt(-2 * log(u1)) * cos(2 * 3.1415926535898 * u2); // Transformation de Box-Muller
    return mu + z * sqrt(sigma);                                   // réalisation du variable aléatoire suivant une loie normale de paramètres mu et sigma
}


double *sfs_infinite(Args args, double **cumulated_branch_lengthes)
{

    double *sfs = calloc(sizeof(double), args.n_samples);
    double size;
    int n_interval = 1 + (args.n_parameters / 2);
    for (int i = 0; i < args.n_samples - 1; i++)
    {
        size = 1.;
        for (int k = 0; k < n_interval; k++)
        {
            sfs[i] += (size * (cumulated_branch_lengthes[i][k + 1] - cumulated_branch_lengthes[i][k]));
            size = args.parameters[args.n_parameters / 2 + k];
        }
        sfs[i] = sfs[i] * args.theta;
        if(args.noised)
            sfs[i] = floor(generate_normal_random(sfs[i], sfs[i])); // sfs.l_genome * sfs.sfs[i] * sfs.theta / 2; // Compute SNPs for each bin
        printf("%f\n", sfs[i]);
    }
    return sfs;
}

<<<<<<< Updated upstream
=======
double *folded_sfs(double *sim_sfs, int sfs_length)
{
    /*
     * Function: folded_sfs
     * ------------------------------------------
     * Folds the simulated SFS data. sfs[i] <- sfs[i] + sfs[n-i] exept if n-i = i
     *
     * Inputs:
     *      sim_sfs : Simulated SFS data.
     *      observed_sfs : Observed SFS data.
     *
     * Output:
     *      folded_sfs : Folded SFS data.
     *
     * Note: This function folds the simulated SFS data to match the observed SFS data by summing corresponding bins.
     */
    double *unfolded = sim_sfs;                                  // Store pointer to unfolded SFS data
    int last_bin = sfs_length / 2;                               // Calculate index of last bin
    int resize = (sfs_length + 1) / 2 + (sfs_length + 1) % 2;    // Calculate size of resized SFS array
    sim_sfs = malloc(sizeof(double) * resize);                   // Allocate memory for resized SFS array
    for (int i = 0; i < (sfs_length + 1) / 2; i++)               // Loop over half of the bins
        sim_sfs[i] = unfolded[i] + unfolded[sfs_length - 1 - i]; // Sum corresponding bins
    if ((sfs_length + 1) % 2 == 0)                               // If number of samples is even
        sim_sfs[last_bin] = sim_sfs[last_bin] / 2.;              // Halve the value of the last bin
    free(unfolded);                                              // Free memory allocated for unfolded SFS data
    return sim_sfs;                                              // Return folded SFS data
}

void write_sfs_to_file(double *sfs, Args args)
{
    //  if(args.oriented)
    //     sfs = folded_sfs(sfs, args.n_samples - 1);
    FILE *file = fopen(args.outfile, "w");
    int replicate = args.noised;
    if (!args.noised)
        replicate = 1;
    if (file == NULL)
    {
        fprintf(stderr, "Error: Unable to open output file.\n");
        exit(1);
    }
    if(args.header)
    {
        fprintf(file, "> sample_size: %d theta: %f parameters:", args.n_samples, args.theta);
        for (int i = 0; i < args.n_parameters; i++)
            fprintf(file, "%f ", args.parameters[i]);
        fprintf(file, "\n");
    }
    for (size_t i = 0; i < replicate; i++)
    {
        for (int i = 0; i < args.n_samples - 1; i++)
            if (args.oriented)
            {
                if (args.noised)
                    fprintf(file, "%.0f ", generate_normal_random(sfs[i], sfs[i])); // sfs.l_genome * sfs.sfs[i] * sfs.theta / 2; // Compute SNPs for each bin
                else
                    fprintf(file, "%.9f ", sfs[i]);
            }
            else
            {
                if (i < args.n_samples / 2)
                {
                    if (args.noised)
                        fprintf(file, "%.0f ", generate_normal_random(sfs[i], sfs[i])); // sfs.l_genome * sfs.sfs[i] * sfs.theta / 2; // Compute SNPs for each bin
                    else
                        fprintf(file, "%.9f ", sfs[i]);
                }
                else
                    fprintf(file, "%d ", 0);
            }
        fprintf(file, "\n ");
    }
    fclose(file);
}

>>>>>>> Stashed changes
int main(int argc, char *argv[])
{
    Args args;
    init_sfs_from_args(argc, argv, &args);

    // Example: Accessing parsed arguments
    printf("n_samples: %d\n", args.n_samples);
    printf("theta: %f\n", args.theta);
    printf("oriented: %d\n", args.oriented);
    printf("seed: %d\n", args.seed);
    printf("outfile: %s\n", args.outfile);
    printf("noised: %d\n", args.noised);

    printf("parameters: ");
    for (int i = 0; i < args.n_parameters; i++)
    {
        printf("%f ", args.parameters[i]);
    }
    printf("\n");

    clock_t start, end;
    double elapsed_time;

    start = clock(); // Début de la mesure du temps
    double *H = absolute_to_relative_times(args.parameters, args.n_parameters);
    double **weight_grid = cumulatve_weight_v2(args.n_samples, args.n_parameters / 2, H);
    double *sfs = sfs_infinite(args, weight_grid);
    end = clock(); // Fin de la mesure du temps
    elapsed_time = ((double)(end - start)) / CLOCKS_PER_SEC; // Conversion en secondes
    printf("Temps d'exécution : %.10f secondes\n", elapsed_time);
    save_cumulated_weight(args.n_samples, args.n_parameters / 2 + 2, weight_grid, "grid.txt");
    // save_cumulated_weight(args.n_samples, args.n_parameters/2 + 1, branch_lengthes, "grid2.txt");
    free_args(&args);
    free(sfs);
    return 0;
}

// 85219043
