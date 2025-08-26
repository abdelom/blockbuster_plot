#ifndef LILNEAR_H
#define LILNEAR_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "blockbuster_grid.h"
#include "sfs.h"

// Definition of the `solution` structure
typedef struct system
{
   double *weight;
   int n_col;
} System; 


// Definition of the `solution` structure
typedef struct solution
{
    double *thetas;           // Dynamic array for population mutation rates between change times, calculated based on the SFS and fixed change times
    int *breakpoints;         // Dynamic array for indices representing times when population size changes occur
    int nb_breakpoints;       // Number of population size changes; 
    double log_likelihood;    // Likelihood of observing the given SFS, based on the scenario with fixed change times and calculated `thetas`
    double distance;          // Distance between the theoretical SFS (simulated with the scenario parameters) and the observed SFS, serving as a goodness-of-fit measure
    double * fitted_sfs;
    double * se_thetas;
    double * residues;
    // double AIC;            // Akaike Information Criterion (AIC), a metric that could aid in model selection by balancing fit and model complexity
} Solution; 



// void regressor_matrix(int n, int m, double *weight, double *regressors);
System init_sytem(int sfs_length, Solution sol, Time_gride tg);
// void replace_negative_with_1(double *theta, size_t n);
// double *SFS_theo(double *thetas, System system, int n)
// double *SFS_to_freq(double *sfs_theo, int n, System system);
// double log_likelihood(double *sfs_obs, System system, int n);
void system_resolution(Solution *sol, SFS sfs, Time_gride tg);
// Solution init_solution();
Solution init_solution_size(int nb_breakpoints);
Solution copy_solution(Solution sol);
void clear_solution(Solution sol);
void save_solution(Solution sol, SFS sfs, Time_gride tg ,char *out_file);
double distance(double *sfs_obs, double *sfs_theo, int n);
double log_likelihood(SFS sfs, System system, Solution *sol, double * sfs_theo);
double *SFS_theo(double *thetas, System system, int n);
void regressor_matrix(int n, int m, double *weight, double *regressors);
void replace_negative_with_1(double *theta, size_t n);
#endif // PAST_FINDER_H
