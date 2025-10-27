#ifndef LILNEAR_H
#define LILNEAR_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "blockbuster_grid.h"
#include "sfs.h"

// Definition of the `System` structure
typedef struct System
{
    double *weight;  // Branch length matrix: each column corresponds to a fixed time interval,
                     // and each row to the number of individuals sustained in the present.
                     // Each entry represents the branch length supporting a given number
                     // of descendants in the present within the associated interval.

    int n_col;       // Number of columns (i.e., number of fixed time intervals),
                     // corresponding to the number of population sizes to estimate
                     // using the least squares estimator.
} System;

// Definition of the `Solution` structure
typedef struct Solution
{
    double *thetas;         // Unknown parameters to estimate: population mutation rates
                            // between change times, obtained using the least squares estimator
                            // based on the observed SFS and fixed change times.

    int *breakpoints;       // Indices representing times when population size changes occur,
                            // delimiting the time intervals where population size is constant.

    int nb_breakpoints;     // Number of population size changes.

    double *time;

    double log_likelihood;  // Log-likelihood of observing the given SFS,
                            // based on the scenario with fixed change times and estimated `thetas`.

    double distance;        // Distance between the expected SFS (computed with the `thetas`
                            // obtained using the least squares estimator) and the observed SFS,
                            // serving as a goodness-of-fit measure.

    double *fitted_sfs;     // Expected SFS computed with the `thetas` obtained using
                            // the least squares estimator.

    double *se_thetas;      // Standard errors of the estimated `thetas`.

    double *residues;       // Residuals between the observed SFS and the expected SFS
                            // computed with the estimated `thetas`.
                            

    // double AIC;           // (Optional) Akaike Information Criterion for model selection.
} Solution;



// void regressor_matrix(int n, int m, double *weight, double *regressors);
System init_sytem(int sfs_length, Solution sol, Time_gride tg);
// void replace_negative_with_1(double *theta, size_t n);
// double *SFS_theo(double *thetas, System system, int n)
// double *SFS_to_freq(double *sfs_theo, int n, System system);
// double log_likelihood(double *sfs_obs, System system, int n);
void system_resolution(Solution *sol, SFS sfs, Time_gride tg);
// Solution init_solution();
Solution init_solution_size(int nb_breakpoints, int gridrefine);
Solution copy_solution(Solution sol);
void clear_solution(Solution sol);
void save_solution(Solution sol, SFS sfs, Time_gride tg , char *out_file, double mut, double gen_time, double genome_length);
double distance(double *sfs_obs, double *sfs_theo, int n);
double log_likelihood(SFS sfs, System system, Solution *sol, double * sfs_theo, Time_gride tg);
double *SFS_theo(double *thetas, System system, int n);
void regressor_matrix(int n, int m, double *weight, double *regressors);
void replace_negative_with_1(double *theta, size_t n);
void convert_times(Solution *sol, Time_gride tg, SFS sfs);
#endif // PAST_FINDER_H
