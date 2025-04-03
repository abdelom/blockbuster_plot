#ifndef PAST_FINDER_H
#define PAST_FINDER_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Definition of the `solution` structure
typedef struct solution
{
    double *thetas;           // Dynamic array for population mutation rates between change times, calculated based on the SFS and fixed change times
    int *breakpoints;         // Dynamic array for indices representing times when population size changes occur
    int nb_breakpoints;       // Number of population size changes; 
    double log_likelihood;    // Likelihood of observing the given SFS, based on the scenario with fixed change times and calculated `thetas`
    double distance;          // Distance between the theoretical SFS (simulated with the scenario parameters) and the observed SFS, serving as a goodness-of-fit measure
    double * se_thetas;
    double * residues;
    double * fitted_sfs;
    // double AIC;            // Akaike Information Criterion (AIC), a metric that could aid in model selection by balancing fit and model complexity
} solution; 


// Fonction pour appliquer une opération de division sur le SFS
double *test_split(double *sfs, int size, int num_blocks);
 
// Fonction pour plier le SFS si non orienté
void fold_sfs(double **sfs, double **cumul_weight, int size, int grid_size);
void sigleton_ignore(double **sfs, double **cumulative_weight, int grid_size);
void singleton_erased(double **sfs, double **cumulative_weight, int sfs_length);

// Fonction pour trouver le scénario
solution* find_scenario(int size, double **cumul_weight, double **sfs, int grid_size, int n_sample, int changes);
solution* recent_infrence(solution *list_solution, int changes, double **sfs, double ** cumul_weight, int sfs_length, int n_sample);
void save_solution(solution sol, int n_sample, char *out_file, double const_ren, int sfs_length);

// Fonction pour libérer la mémoire allouée pour l'intégrale des poids
void free_integral_grid(double **cumul_weight, int n_sample);

// Fonction pour libérer la mémoire allouée pour la solution
void clear_solution(solution sol);

#endif // PAST_FINDER_H
