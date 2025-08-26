#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <lapacke.h>
#include <stddef.h>
#include <cblas.h>
#include <omp.h>
#include "sfs.h"
#include "linear_regression.h"
#include "linear_regression_f.h"
#include "blockbuster_grid.h"



Flag init_flag(int sfs_length, double *thetas , double n_theta)
{
    Flag flag;
    flag.thetas = thetas;
    flag.n_theta = n_theta;
    flag.sfs = malloc(sizeof(double) * sfs_length);
    flag.n_theta_fixed = 0;
    for(int i = 0; i < n_theta; i ++)
    {
        if(thetas[i] != 0)
            flag.n_theta_fixed ++;
    }
    return flag;
}


double * init_system_f(SFS sfs, Flag *flag, System system)
{
    // // System system = init_sytem(sfs.sfs_length, sol, tg);
    double * w = calloc((system.n_col - flag->n_theta_fixed) * sfs.sfs_length, sizeof(double));
    for(int k = 0; k < sfs.sfs_length; k ++)
        flag->sfs[k] = sfs.training[k];
    int jw = 0;
    for(int j = 0; j < system.n_col; j ++)
    {
        if(flag->thetas[j] > 0.)
        {
            for(int k = 0; k < sfs.sfs_length; k ++)
                flag->sfs[k] -= flag->thetas[j] * system.weight[k * system.n_col + j];
                // system.weight[k * system.n_col + j] = -1.;
        }
        else
        {
            for(int k = 0; k < sfs.sfs_length; k ++)
                w[k * (system.n_col - flag->n_theta_fixed) + jw] = system.weight[k * system.n_col + j];
            jw += 1;
        }
    }
    // int nb_0 = 0;
    // for(int i = 0; i + nb_0 < sfs.sfs_length * system.n_col; i ++)
    // {
    //     while(system.weight[i + nb_0] < 0.)
    //         nb_0 += 1;
    //     system.weight[i] = system.weight[i + nb_0];
    // }
    // system.weight = realloc(system.weight, sfs.sfs_length * system.n_col - nb_0);
    // system.n_col -= nb_0 / sfs.sfs_length;
    return w;
}


/**
 * Resolves the system of equations to estimate population mutation rates (thetas) and calculates log likelihood and distance from observed sfs given fixed times of change.
 *
 * @param sol         Pointer to the solution structure containing parameters for the regression.
 * @param sfs         Pointer to a 2D array representing the Site Frequency Spectra (SFS).
 * @param cumul_weight Pointer to a 2D array representing cumulative weights for regression.
 * @param sfs_length  Length of the SFS, representing the number of observations.
 *
 * This function performs the following tasks:
 * Assembles regression weights based on the cumulative branch lengths.
 * Performs least square method to estimate the population mutation rates (thetas).
 **/
void system_resolution_f(Solution *sol, SFS sfs, Time_gride tg, Flag flag)
{
    // Step 1: Assemble regression weights from cumulative branch lengths
    System system = init_sytem(sfs.sfs_length, *sol, tg);
    // double *weight_f = init_system_f(sfs, &flag, system);
    // double *regressors = calloc(system.n_col * sfs.sfs_length, sizeof(double));
    int n_col = flag.n_theta - flag.n_theta_fixed;
    double *thetas_tmp = calloc(n_col, sizeof(double));
    // Step 3: Compute the regression matrix using the regression weights (X^T X)-1X^T
    if(flag.n_theta > flag.n_theta_fixed)
    {
        double *weight_f = init_system_f(sfs, &flag, system);
        double *regressors = calloc(n_col * sfs.sfs_length, sizeof(double));
        regressor_matrix(sfs.sfs_length, n_col, weight_f, regressors);
        // Step 4: Estimate population mutation rates (thetas) using matrix multiplication (X^T X)-1X^T * SFS
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    n_col, 1, sfs.sfs_length, 1.0, regressors, sfs.sfs_length, flag.sfs, 1, 0.0, thetas_tmp, 1);
        replace_negative_with_1(thetas_tmp, n_col); // if thetas are negatives they are replaced by 1 as population sizes cannot be inferior to 0}
        free(regressors);
        free(weight_f);
    }
    int i_f = 0, i_th = 0;
    for (int i = 0; i < system.n_col; i++)
    {
        // printf("%f %f %d %d %d\n", flag.thetas[i_f] , thetas_tmp[0], i_f, i_th, (flag.thetas[i_f] > 0.));
        sol->thetas[i] = (flag.thetas[i_f] > 0.) ? flag.thetas[i_f] : thetas_tmp[i_th++];
        i_f +=1;
        // printf("%f\n", sol->thetas[i]);
    }
    double *sfs_theo = SFS_theo(sol->thetas, system, sfs.sfs_length);
    sol->log_likelihood = log_likelihood(sfs, system, sol, sfs_theo);
    sol->distance = distance(sfs.training, sfs_theo, sfs.sfs_length);
    free(system.weight);
    free(sfs_theo);    // Free the memory allocated for the frequency SFS
    free(thetas_tmp);
}


