#include "sfs.h"
#include "blockbuster_grid.h"
#include "linear_regression.h"
#include "blockbuster.h"
#include <stdlib.h>  // pour malloc, free, qsort

// Prototype du comparateur pour qsort
int compare_double(const void *a, const void *b) {
    double diff = *(const double *)a - *(const double *)b;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    return 0;
}

double *adapt_scale(double * time_scale, Solution sol, int repeats)
{
    double * thetas = malloc((repeats * sol.nb_breakpoints + 1) * sizeof(double));
    int j = 0;
    for(int i = 0; i < repeats * sol.nb_breakpoints + 1; i ++)
    {
        if (j < sol.nb_breakpoints)
        { 
            if(sol.time[j] == time_scale[i])
                j ++;
        }
        thetas[i] = sol.thetas[j];
    }
    return thetas;
}

// Fonction pour calculer la moyenne harmonique par ligne
double *harmonic_mean_rows(double **matrix, int n_rows, int n_cols)
{
    double *result = malloc(n_cols * sizeof(double));
    for (int i = 0; i < n_cols; i++) {
        double sum_inverse = 0.0;
        int count = 0;

        for (int j = 0; j < n_rows; j++) {
            // if (matrix[i][j] != 0.0) {  // éviter la division par zéro
                sum_inverse += 1.0 / matrix[j][i ];
                count++;
            // }
        }

        // if (count > 0)
            result[i] = count / sum_inverse;
        // else
        //     result[i] = 1.0;  // ou NAN si tu préfères
    }
    return result;
}

int format_bootstrap(Solution *list_solution, int repeats, double mutation_rate, double genome_length, double generation_time, const char *prefix)
{
    int nb_solutions = repeats;
    int nb_breakpoints = list_solution[0].nb_breakpoints;
    if(repeats < 20)
        return 0;
    // Allocation du tableau
    double *time_scale = malloc(sizeof(double) * nb_solutions * nb_breakpoints);
    double **rep = malloc(sizeof(double*) * nb_solutions);
    if (!time_scale || !rep) {
        fprintf(stderr, "Erreur d’allocation mémoire pour time_scale ou rep\n");
        exit(EXIT_FAILURE);
    }

    // Copie des valeurs de temps dans le tableau
    for (int i = 0; i < nb_solutions; i++) {
        for (int j = 0; j < nb_breakpoints; j++) {
            time_scale[i * nb_breakpoints + j] = list_solution[i].time[j];
        }
    }

    // Tri du tableau par ordre croissant
    qsort(time_scale, nb_solutions * nb_breakpoints, sizeof(double), compare_double);

    // Adapter chaque solution à la nouvelle échelle
    for (int i = 0; i < nb_solutions; i++) {
        rep[i] = adapt_scale(time_scale, list_solution[i], repeats);
    }

    // Calcul de la moyenne harmonique
    double *har_m = harmonic_mean_rows(rep, repeats, repeats * list_solution[0].nb_breakpoints + 1);

    // Construction du chemin du fichier de sortie
    char filepath[512];
    snprintf(filepath, sizeof(filepath), "%s/boot.txt", prefix);

    FILE *f = fopen(filepath, "w");
    if (!f) {
        fprintf(stderr, "Erreur : impossible de créer le fichier %s\n", filepath);
        free(time_scale);
        free(rep);
        free(har_m);
        exit(EXIT_FAILURE);
    }

    // Écriture dans le fichier : time_scale[i] et har_m[i]
    int total = nb_solutions * nb_breakpoints;
    fprintf(f, "time_in_Ne0_generation \t Theta");
    if(mutation_rate > 0 && genome_length > 0)
    {
        fprintf(f, "\t time_in_generation \t Effective size");
        if(generation_time > 0)
            fprintf(f, "\t time_in_years");
    }
    for (int i = 0; i < total; i++)
    {
        fprintf(f, "\n");
        fprintf(f, "%f \t %f ", time_scale[i], har_m[i]);
        if(mutation_rate > 0 && genome_length > 0)
        {
            fprintf(f, "\t %f \t %f ", 2 * time_scale[i] * har_m[0] / (4 * mutation_rate * genome_length), har_m[i] / (4 * mutation_rate * genome_length));
            if(generation_time > 0)
                 fprintf(f, "\t %f ", 2 * generation_time * time_scale[i] * har_m[0] / (4 * mutation_rate * genome_length));

        }
    }
    fclose(f);

    // Libération de la mémoire
    free(time_scale);
    for (int i = 0; i < nb_solutions; i++) {
        free(rep[i]);
    }
    free(rep);
    free(har_m);
    for(int b = 0; b < repeats; b ++)
    {
        free(list_solution[b].time);
        clear_solution(list_solution[b]);
    }
    return 1;
}


// comparaison pour qsort
int compare_ints(const void *a, const void *b) {
    int x = *(const int *)a;
    int y = *(const int *)b;
    return (x > y) - (x < y);
}

int* tirer_uniques_tries(int n, int k) {
    if (n > k) {
        // impossible de tirer n entiers distincts entre 1 et k
        return NULL;
    }

    int *tab = malloc(n * sizeof(int));
    if (!tab) return NULL;

    int count = 0;

    while (count < n) {
        int val = rand() % k + 1;

        // vérifier si val existe déjà
        int existe = 0;
        for (int i = 0; i < count; i++) {
            if (tab[i] == val) {
                existe = 1;
                break;
            }
        }

        if (!existe) {
            tab[count++] = val;
        }
    }

    // tri final
    qsort(tab, n, sizeof(int), compare_ints);

    return tab;
}


int recursive_bk_combination_sampling(int *brk, int nb_breakpoints, int grid_size)
{
    // Start with the last breakpoint in the array
    int i = nb_breakpoints - 1;
    int *comb = malloc(nb_breakpoints * sizeof(int));
    // Attempt to increment the last breakpoint by 1
    brk[i] += 1;

    // Adjust previous breakpoints if the current one exceeds allowable range
    while (i >= 0 && brk[i] > grid_size + 1 - nb_breakpoints + i)
    {
        // Move to the previous breakpoint in the array
        i--;
        // Increment the previous breakpoint
        brk[i] += 1;
    }
    // For any indices after `i`, reset them to maintain increasing order
    for (int j = i + 1; j < nb_breakpoints; j++)
        brk[j] = brk[j - 1] + 1;
    // If the first breakpoint reaches its upper limit, all combinations have been generated

    if (i == 0 && brk[0] == grid_size + 1 - nb_breakpoints)
        return 0;
    // Return 1 to indicate more combinations are available
    return 1;
}

// Solution generate_brk_combinations_sampling(int nb_breakpoints, SFS sfs, Time_gride tg)
// {
//     // Initialize the solution with specified number of breakpoints
//     Solution sol = init_solution_size(nb_breakpoints, 1);

//     // Set the last breakpoint to the grid limit (end boundary)
//     sol.breakpoints[nb_breakpoints] = GRIDREFINE * tg.grid_size + 2;
//     int * sampled_brk = tirer_uniques_tries(tg.grid_size, GRIDREFINE * tg.grid_size + 1);
//     int * current_combination = malloc(nb_breakpoints * sizeof(int));
//     for(int i = 0; i < nb_breakpoints; i ++)
//     {
//         current_combination[i] = sol.breakpoints[i];
//         sol.breakpoints[i] = sampled_brk[current_combination[i]];
//     }
//     // Flag for stopping the combination generation
//     int arret = 1;
//     // Calculate the initial solution with the given breakpoints
//     system_resolution(&sol, sfs, tg);
//     // Initialize a temporary solution to keep track of the best found solution
//     Solution tmp_sol = copy_solution(sol);
//     // Loop through all breakpoint combinations until `arret` is set to 0
//     while (arret && nb_breakpoints > 0)
//     {
//         // Generate the next combination of breakpoints
//         arret = recursive_bk_combination_sampling(current_combination, sol.nb_breakpoints, tg.grid_size);
//         for(int i = 0; i < nb_breakpoints; i ++)
//             sol.breakpoints[i] = sampled_brk[current_combination[i] - 1];
//         // Calculate the solution (log-likelihood and distance) for the new combination
//         system_resolution(&sol, sfs, tg);
//         // If the new solution has a higher log-likelihood and all theta values are positive, keep its
//         if (isnan(tmp_sol.log_likelihood) || sol.log_likelihood > tmp_sol.log_likelihood)
//         {
//             clear_solution(tmp_sol);
//             tmp_sol = copy_solution(sol);
//         }
//     }
//     clear_solution(sol);
//     return tmp_sol;
// }

Solution * bootstrap(SFS sfs, Time_gride tg, int epochs, int repeats)
{
    Solution *list_solution = malloc(repeats  * sizeof(Solution));
    // SFS sfs_boot = noise_sfs(sfs);
    for (int b = 0; b < repeats; b++)
    {
        // printf("\n x x ");
        SFS sfs_boot = noise_sfs(sfs);
        // Solution *solution = find_scenario(sfs_boot, tg, epochs);
        sfs_boot.delta_time = 0.1;       // char *outfile_b = construct_output_filepath(out_file, "bootstrap_scenarios.txt");
        list_solution[b] = generate_brk_combinations(epochs - 1, sfs_boot, tg);
        if (epochs - 1 >= 1)
        {
            refine_solution(&list_solution[b], sfs_boot, tg);
        }
        if((b + 1) % 50 == 0)
            printf("\n smoothing completion : %f % ", 100 * (float)(b + 1) / (float)repeats);
        convert_times(&list_solution[b], tg, sfs_boot);
        // save_list_solution(list_solution, sfs_boot, outfile_b, epochs, tg);
        clear_sfs(sfs_boot);
        // clear_solution(sol);
    }
    // format_bootstrap(list_solution);
    return list_solution;
}