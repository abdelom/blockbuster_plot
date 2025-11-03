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
#include "blockbuster_grid.h"


/**
 * Computes regression weights representing branch lengths for lineages with a certain number of descendants
 * in the present, within specified time intervals, based on cumulative branch lengths
 * in `integral_grid`.
 *
 * @param sfs_length   Length of the Site Frequency Spectrum (SFS).
 * @param n_theta      Number of fixed time intervals (and corresponding mutation rates).
 * @param time_index   Array of indices marking upper bounds of each fixed time interval.
 * @param integral_grid 2D array where each entry `integral_grid[k][i]` represents the cumulative branch
 *                      length for `k + 1` descendants in the present, through time up to index `i`.
 *
 * @return A dynamically allocated array `C_kj_matrix` of size `n_theta * sfs_length` containing
 *         the branch lengths (weights) for each time interval and sample frequency. The caller is
 *         responsible for freeing this memory.
 *
 * Example usage:
 *     double *weights = weight_assembly_1d(sfs_length, n_theta, time_index, integral_grid);
 *     // Use `weights` array here
 *     free(weights);  // Free the allocated memory when done
 */
System init_sytem(int sfs_length, Solution sol, Time_gride tg)
{
    System system;
   
    system.n_col =  (sol.nb_breakpoints + 1);
    system.weight = calloc(system.n_col * sfs_length, sizeof(double));
    // Iterate over each number of descendants in the present, up to the length of the SFS
    for (int k = 0; k < sfs_length; k++)
    {
        int lower_time_index = 0; // Initialize the lower bound of the time interval to 0

        // For each time interval, compute the weight based on cumulative branch length differences
        for (int j = 0; j < system.n_col; j++)
        {
            int upper_time_index = sol.breakpoints[j]; // Upper bound of the current time interval

            // Calculate the weight by taking the difference in cumulative branch lengths, representing
            // the branch length for `k + 1` descendants within the current time interval
            system.weight[k * system.n_col + j] = tg.cumulative_bl[k][upper_time_index] - tg.cumulative_bl[k][lower_time_index];

            // Move the lower bound to the current upper bound for the next interval
            lower_time_index = upper_time_index;
        }
    }
    return system; // Return the computed regression weight matrix
}


double varcovar(int n, int m, double *weight, double *XXT)
{
    // Compute the product X * X^T
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                m, m, n, 1.0, weight, m, weight, m, 0.0, XXT, m);
    lapack_int ipiv[m]; // Array to store pivot indices for LU decomposition
                        // Perform LU decomposition on the Gram matrix XXT
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, m, XXT, m, ipiv);
    // Invert the Gram matrix using the LU decomposition
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, m, XXT, m, ipiv);
}


// Function to perform linear regression: (X * X^T) * X^T * Y
/**
 * Function to perform linear regression using the formula: (X * X^T) * X^T * Y
 *
 * @param n         Length of the Site Frequency Spectrum (SFS), representing the number of observations.
 * @param m         Number of population size changes (features).
 * @param weight    Pointer to the array representing the regression weights (branch lengths for a certain number of descendants in the present).
 * @param regressors Pointer to the array where the computed regression coefficients (result) will be stored.
 *
 * This function computes the regression coefficients using the Ordinary Least Squares method:
 * 1. It computes the matrix product X * X^T to obtain the Gram matrix.
 * 2. It performs LU decomposition of the Gram matrix to prepare for inversion.
 * 3. It inverts the Gram matrix to solve for the regression coefficients.
 * 4. Finally, it multiplies the inverted Gram matrix with X^T and Y to obtain the regression coefficients.
 */
void regressor_matrix(int n, int m, double *weight, double *regressors)
{
    double *XXT = calloc(m * m, sizeof(double)); // Array to store the Gram matrix (X * X^T)
    // Compute the product X * X^T
    varcovar(n, m, weight, XXT);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                m, n, m, 1.0, XXT, m, weight, m, 0.0, regressors, n);
    free(XXT);
}

// Replaces all negative values in a population mutation rate vector with 1.
//
// Parameters:
// - theta: Pointer to a vector of population mutation rates obtained via 
//          a least-squares method. The parameter space is not constrained, 
//          so negative values may occur.
// - n: Size of the vector.
//
// This function iterates through the vector `theta`. If any element is 
// negative, it is replaced with 1. This ensures that the vector contains 
// only valid mutation rates.
//
// Arguments:
// - double *theta: A vector of population mutation rates.
// - size_t n: The number of elements in the vector.
//
// Complexity: O(n), where n is the size of the vector.
// void replace_negative_with_1(double *theta, size_t n)
// {
//     for (size_t i = 0; i < n; i++) // Iterate through each element of the vector.
//     {
//         if (theta[i] < 1e4) // Check if the current value is negative.
//             theta[i] = 5e8; // Replace the negative value with 1.
//     }
// }



/**
 * Computes the theoretical Site Frequency Spectrum (SFS) based on estimated population mutation rates and regression weights.
 *
 * @param thetas    Pointer to the array of population mutation rates for each time interval.
 * @param weight    Pointer to the array representing the regression weights (branch lengths for a certain number of descendants in the present).
 * @param n         Length of the Site Frequency Spectrum (SFS), representing the number of observations.
 * @param m         Number of population size changes (features).
 *
 * @return A dynamically allocated array `sfs_theo` of size `n` containing the computed theoretical SFS.
 *         The caller is responsible for freeing this memory.
 *
 * This function uses matrix multiplication to compute the theoretical SFS:
 * It performs the operation SFS_theo = weight * thetas, where:
 * - `weight` is an m x n matrix (number of changes x observations),
 * - `thetas` is an m x 1 matrix (number of changes x 1),
 * - The result is an n x 1 matrix representing the theoretical SFS.
 */
double *SFS_theo(double *thetas, System system, int n)
{
    // Allocate memory for the theoretical SFS array
    double *sfs_theo = malloc(sizeof(double) * n);
    if (sfs_theo == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Perform matrix multiplication to compute the theoretical SFS: sfs_theo = weight * thetas
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, system.n_col, 1.0, system.weight, system.n_col, thetas, 1, 0.0, sfs_theo, 1);

    return sfs_theo; // Return the computed theoretical SFS
}

/**
 * Converts a Site Frequency Spectrum (SFS) from counts to frequencies.
 *
 * @param sfs_theo Pointer to the array representing the theoretical SFS (counts).
 * @param n        Length of the SFS, representing the number of observations.
 *
 * @return A dynamically allocated array `sfs_theo2` of size `n` containing the computed frequency SFS.
 *         The caller is responsible for freeing this memory.
 *
 * This function computes the frequency of each site in the SFS by dividing the counts
 * by the total number of sites (n_snp). The frequencies are calculated as follows:
 */
double *SFS_to_freq(double *sfs_theo, System system, int sfs_length)
{
    double n_snp = 0;                               // Variable to hold the total number of SNPs
    double *sfs_theo2 = malloc(sizeof(double) * sfs_length); // Allocate memory for the frequency array
    // Calculate the total number of SNPs by summing the counts in sfs_theo
    for (int i = 0; i < sfs_length; i++)
        n_snp += sfs_theo[i];
    // Compute the frequency for each site and store in sfs_theo2
    for (int i = 0; i < sfs_length; i++)
        sfs_theo2[i] = sfs_theo[i] / n_snp;
    return sfs_theo2; // Return the computed frequency SFS
}

/**
 * Calculates the log likelihood of the observed Site Frequency Spectrum (SFS) given the theoretical SFS, computed from etimated parameters.
 *
 * @param sfs_obs Pointer to the array representing the observed SFS (counts).
 * @param sfs_theo Pointer to the array representing the theoretical SFS (counts).
 * @param n        Length of the SFS, representing the number of observations.
 *
 * @return The log likelihood value computed from the observed and theoretical SFS.
 *
 * This function computes the log likelihood using the formula:
 * L = sum(sfs_obs[i] * log(sfs_theo2[i])), where sfs_theo2 is the frequency representation of sfs_theo.
 * The process involves the following steps:
 * 1. Convert the theoretical SFS from counts to frequencies using the SFS_to_freq function.
 * 2. Compute the log likelihood by summing the product of the observed counts and the logarithm of the theoretical frequencies.
 */
double log_likelihood(SFS sfs, System system, Solution *sol, double * sfs_theo, Time_gride tg)
{
    double llikelihood = 0;                       // Variable to hold the calculated log likelihood
    double *sfs_theo2 = SFS_to_freq(sfs_theo, system, sfs.sfs_length); // Convert theoretical SFS to frequencies
    // Calculate the log likelihood
    double time = 0.0, lb = 0.0, delta_time = 0.0;
    if(sol->thetas[sol -> nb_breakpoints] < 1e4)
        return -INFINITY;
    if(sol->thetas[0] < 1e4)
        return -INFINITY;
    for (int i = 0; i < sol -> nb_breakpoints; i ++)
    {
        delta_time = (tg.time_scale[sol->breakpoints[i] - 1] - lb) * sol->thetas[i] / sol->thetas[0];
        // if(delta_time < 0)
        
        if(time >= 0. && (delta_time / time < sfs.delta_time)) 
        {
            return -INFINITY;
        }
        time += delta_time ;
        lb = tg.time_scale[sol->breakpoints[i] - 1];

    }
    for (int i = 0; i < sfs.sfs_length; i++)
        llikelihood += (sfs.test[i] * log10(sfs_theo2[i]));
    // for (int i = 1; i <= sol -> nb_breakpoints; i ++)
        // llikelihood += - 1 * fabs(log10(sol->thetas[i] / sol->thetas[i-1]));
    free(sfs_theo2);
    return llikelihood; // Return the computed log likelihood
}


/**
 * Computes the Euclidean distance between the observed and theoretical Site Frequency Spectra (SFS).
 *
 * @param sfs_obs Pointer to the array representing the observed SFS (counts).
 * @param sfs_theo Pointer to the array representing the theoretical SFS (counts).
 * @param n        Length of the SFS, representing the number of observations.
 *
 * @return The computed Euclidean distance between the observed and theoretical SFS.
 */
double distance(double *sfs_obs, double *sfs_theo, int n)
{
    double distance = 0; // Variable to hold the accumulated squared differences

    // Calculate the sum of squared differences between observed and theoretical SFS
    for (int i = 0; i < n; i++)
    {
        distance += powf((sfs_obs[i] - sfs_theo[i]), 2.); // Squared difference
    }
    return sqrt(distance); // Return the square root of the accumulated distance
}


/**
 * Calculates the standard errors (confidence intervals) for the estimated 
 * regression weights (thetas).
 *
 * @param sol          Pointer to the solution structure containing thetas and other parameters.
 * @param sfs_length   Length of the Site Frequency Spectrum (number of observations).
 * @param cumul_weight Pointer to a 2D array of cumulative weights for regression.
 *
 * Steps:
 * 1. Compute the variance of errors (error_var) based on the distance and degrees of freedom.
 * 2. Assemble regression weights and compute the variance-covariance matrix.
 * 3. Extract diagonal elements of the variance-covariance matrix to calculate standard errors.
 */
void thetas_se(Solution *sol, int sfs_length, Time_gride tg)
{
    int nb_thetas = sol->nb_breakpoints + 1;
    // Calculate error variance using the residual sum of squares
    double error_var = sol->distance * sol->distance / (sfs_length - nb_thetas);
    // Assemble regression weights
    System system = init_sytem(sfs_length, *sol, tg);
    // double *weight = weight_assembly_1d(sfs_length, nb_thetas, sol->breakpoints, cumul_weight);

    // Allocate memory for variance-covariance matrix
    double *XXT = calloc(nb_thetas * nb_thetas, sizeof(double));
    sol->se_thetas = calloc(nb_thetas, sizeof(double));

    // Compute variance-covariance matrix
    varcovar(sfs_length, nb_thetas, system.weight, XXT);

    // Calculate standard errors for each theta
    for (int j = 0; j < nb_thetas; j++)
        sol->se_thetas[j] = sqrt(error_var * XXT[j * nb_thetas + j]);

    // Free allocated memory
    free(system.weight);
    free(XXT);
}

// Calculates the residuals of the linear regression for the observed SFS.
//
// Parameters:
// - sol: Pointer to the solution structure containing breakpoints and thetas.
// - cumul_weight: 2D array of cumulative weights for the regression.
// - sfs_length: Length of the observed SFS.
// - sfs: Array of observed SFS values.
//
// Steps:
// 1. Compute regression weights using `weight_assembly_1d`.
// 2. Calculate the theoretical SFS (`sol->fitted_sfs`) from the weights and thetas.
// 3. Compute residuals as the difference between observed and theoretical SFS.
void residues(Solution *sol, Time_gride tg, SFS sfs)
{
    int nb_thetas = sol->nb_breakpoints + 1;
    System system = init_sytem(sfs.sfs_length, *sol, tg);
    sol->fitted_sfs = SFS_theo(sol->thetas, system, sfs.sfs_length);
    sol->residues = calloc(sfs.sfs_length, sizeof(double));
    for (int i = 0; i < sfs.sfs_length; i++)
        sol->residues[i] = sfs.training[i] - sol->fitted_sfs[i];
    free(system.weight);
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
void system_resolution(Solution *sol, SFS sfs, Time_gride tg)
{
    // Step 1: Assemble regression weights from cumulative branch lengths
    System system = init_sytem(sfs.sfs_length, *sol, tg);
    double *regressors = calloc(system.n_col * sfs.sfs_length, sizeof(double));
    // Step 3: Compute the regression matrix using the regression weights (X^T X)-1X^T
    regressor_matrix(sfs.sfs_length, system.n_col, system.weight, regressors);
    // Step 4: Estimate population mutation rates (thetas) using matrix multiplication (X^T X)-1X^T * SFS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                system.n_col, 1, sfs.sfs_length, 1.0, regressors, sfs.sfs_length, sfs.training, 1, 0.0, sol->thetas, 1);
    // replace_negative_with_1(sol->thetas, system.n_col); // if thetas are negatives they are replaced by 1 as population sizes cannot be inferior to 0
    double *sfs_theo = SFS_theo(sol->thetas, system, sfs.sfs_length);
    sol->log_likelihood = log_likelihood(sfs, system, sol, sfs_theo, tg);
    sol->distance = distance(sfs.training, sfs_theo, sfs.sfs_length);
    free(system.weight);
    free(sfs_theo);    // Free the memory allocated for the frequency SFS
    free(regressors);
}


// Function to initialize a solution with a specified number of breakpoints
Solution init_solution_size(int nb_breakpoints, int gridrefine)
{
    Solution sol;
    sol.nb_breakpoints = nb_breakpoints;
    // Allocate memory for `nb_breakpoints + 1` to include initial and final states or boundaries, the last one corespond to the infinity
    sol.breakpoints = calloc(sizeof(int), (nb_breakpoints + 1));
    sol.thetas = calloc(sizeof(double), (nb_breakpoints + 1)); // Allocate space for each interval’s mutation rate

    for (int i = 0; i < nb_breakpoints; i++){
        sol.breakpoints[i] = i * gridrefine + 1; // Initialize breakpoints as a simple sequence; this may be customized to represent actual change times
    }
        
    return sol;
}


// Function to copy an existing solution into a new structure
Solution copy_solution(Solution sol)
{
    Solution solution_c;
    solution_c.log_likelihood = sol.log_likelihood;                          // Copy the likelihood from the original solution
    solution_c.distance = sol.distance;                                      // Copy the distance from the original solution
    solution_c.nb_breakpoints = sol.nb_breakpoints;                          // Copy the number of breakpoints
    solution_c.breakpoints = malloc(sizeof(int) * (sol.nb_breakpoints + 1)); // Allocate space for breakpoints in the new copy
    solution_c.thetas = malloc(sizeof(double) * (sol.nb_breakpoints + 1));   // Allocate space for mutation rates in the new copy
    for (int i = 0; i < sol.nb_breakpoints + 1; i++)
    {
        solution_c.breakpoints[i] = sol.breakpoints[i]; // Copy each breakpoint index from the original
        solution_c.thetas[i] = sol.thetas[i];           // Copy each mutation rate from the original
    }
    solution_c.se_thetas = NULL;
    return solution_c;
}

// Function to free the memory associated with a solution
void clear_solution(Solution sol)
{
    free(sol.thetas); // Free memory allocated for mutation rates
    // if (sol.breakpoints)
    free(sol.breakpoints); // Free memory allocated for breakpoints, if they were allocated
}


/**
 * Saves the solution if all mutation rates are positive.
 *
 * Outputs solution details, including log-likelihood, distance between theoretical and observed SFS,
 * time indices (breakpoints), mutation rates (thetas), and absolute times at each breakpoint. The
 * absolute times are derived by first calculating relative times from indices, then adjusting based on
 * effective population size changes.
 *
 * @param sol      The solution structure containing model parameters.
 * @param n_sample Sample size, used for normalizing the time values.
 */
void save_solution(Solution sol, SFS sfs, Time_gride tg , char *out_file, double mut, double gen_time, double genome_length)
{
    // Ensure the directory exists before proceeding
    double * effective_Ne = calloc(sol.nb_breakpoints + 1, sizeof(double));
    FILE *file;
    if(sol.nb_breakpoints > 0)
        file = fopen(out_file, "a");
    else
        file = fopen(out_file, "w");
    thetas_se(&sol, sfs.sfs_length, tg);
    residues(&sol, tg, sfs);
    fprintf(file, " > %d epochs model \n log_likelihood: %f ", sol.nb_breakpoints + 1, sol.log_likelihood);
    fprintf(file, "\n distance: %f ", sol.distance);
    fprintf(file, "\n thetas: ");
    for (int i = 0; i <= sol.nb_breakpoints; i++)
        fprintf(file, "%f ", sol.thetas[i] / sfs.training_size);
    if(mut > 0 && genome_length > 0){
        fprintf(file, "\n Effective size: ");
        for(int i = 0; i <= sol.nb_breakpoints; i++){
            effective_Ne[i] = sol.thetas[i] / (4. * mut * genome_length * sfs.training_size);
            fprintf(file, "%f ", effective_Ne[i]);
        }
    }
    fprintf(file, "\n time in unit of Ne generations: ");
    for (int i = 0; i < sol.nb_breakpoints; i++)
        fprintf(file, "%f ", sol.time[i]);
    if(mut > 0 && genome_length > 0){
        fprintf(file, "\n time in generations: ");
        for(int i = 0; i < sol.nb_breakpoints; i++){
             fprintf(file, "%f ", sol.time[i] * effective_Ne[0] * 2);
        }
         if(gen_time > 0){
        fprintf(file, "\n time in years: ");
        for(int i = 0; i < sol.nb_breakpoints; i++){
             fprintf(file, "%f ", sol.time[i] * effective_Ne[0] * gen_time * 2);
        }
    }
    }
    if (sol.se_thetas != NULL)
    {
        fprintf(file, "\n se_thetas: ");
        for (int i = 0; i <= sol.nb_breakpoints; i++)
            fprintf(file, "%f ", sol.se_thetas[i] / sfs.training_size);
    }
    if (sol.residues != NULL)
    {
        fprintf(file, "\n residues: ");
        for (int i = 0; i < sfs.sfs_length; i++)
            fprintf(file, "%f ", sol.residues[i] / sfs.training_size);
    }
    if (sol.fitted_sfs != NULL)
    {
        fprintf(file, "\n fitted_sfs: ");
        for (int i = 0; i < sfs.sfs_length; i++)
            fprintf(file, "%f ", sol.fitted_sfs[i] / sfs.training_size);
    }
    fprintf(file, "\n");
    fclose(file);
    free(sol.residues);
    free(sol.se_thetas);
    free(sol.fitted_sfs);
    free(sol.time);    // clear_solution(sol);
}

void convert_times(Solution *sol, Time_gride tg, SFS sfs)
{
    double time = 0.0, lb = 0.0, delta_time = 0.0;
    sol-> time = calloc(sol->nb_breakpoints, sizeof(double));
    for (int i = 0; i < sol -> nb_breakpoints; i ++)
    {
        delta_time = (tg.time_scale[sol->breakpoints[i] - 1] - lb) * sol->thetas[i] / sol->thetas[0];
        time += delta_time ;
        lb = tg.time_scale[sol->breakpoints[i] - 1];
        sol -> time[i] = time;
    }
}