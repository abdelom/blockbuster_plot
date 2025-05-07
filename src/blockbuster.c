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
#include "blockbuster.h"

// Function to initialize an empty solution with minimal allocation
solution init_solution()
{
    solution sol;
    sol.thetas = malloc(sizeof(double)); // Allocate memory for a single mutation rate (could be resized later)
    sol.breakpoints = NULL;              // No breakpoints initially, as no changes are assumed
    sol.nb_breakpoints = 0;              // Start with zero size changes, making this an empty or initial solution
    sol.log_likelihood = 0;              // Initialize likelihood to zero, to be calculated after SFS fitting
    sol.distance = 0;                    // Initialize distance to zero, to be calculated later
    sol.se_thetas = NULL;
    sol.residues = NULL;
    sol.fitted_sfs = NULL;
    // sol.AIC = 0;                           // Optional: AIC initialization, could be implemented for model comparison
    return sol;
}

// Function to initialize a solution with a specified number of breakpoints
solution init_solution_size(int nb_breakpoints)
{
    solution sol;
    sol.nb_breakpoints = nb_breakpoints;
    // Allocate memory for `nb_breakpoints + 1` to include initial and final states or boundaries, the last one corespond to the infinity
    sol.breakpoints = calloc(sizeof(int), (nb_breakpoints + 1));
    sol.thetas = calloc(sizeof(double), (nb_breakpoints + 1)); // Allocate space for each interval’s mutation rate
    sol.se_thetas = NULL;
    sol.residues = NULL;
    for (int i = 0; i < nb_breakpoints; i++)
        sol.breakpoints[i] = i * GRIDREFINE + 1; // Initialize breakpoints as a simple sequence; this may be customized to represent actual change times
    return sol;
}

// Function to copy an existing solution into a new structure
solution copy_solution(solution sol)
{
    solution solution_c;
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
void clear_solution(solution sol)
{
    free(sol.thetas); // Free memory allocated for mutation rates
    // if (sol.breakpoints)
    free(sol.breakpoints); // Free memory allocated for breakpoints, if they were allocated
}

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
double *weight_assembly_1d(int sfs_length, int n_theta, int *time_index, double **integral_grid)
{
    // Allocate memory for the regression weight matrix, initializing values to 0
    double *C_kj_matrix = calloc(n_theta * sfs_length, sizeof(double));

    // Iterate over each number of descendants in the present, up to the length of the SFS
    for (int k = 0; k < sfs_length; k++)
    {
        int lower_time_index = 0; // Initialize the lower bound of the time interval to 0

        // For each time interval, compute the weight based on cumulative branch length differences
        for (int j = 0; j < n_theta; j++)
        {
            int upper_time_index = time_index[j]; // Upper bound of the current time interval

            // Calculate the weight by taking the difference in cumulative branch lengths, representing
            // the branch length for `k + 1` descendants within the current time interval
            C_kj_matrix[k * n_theta + j] = integral_grid[k][upper_time_index] - integral_grid[k][lower_time_index];

            // Move the lower bound to the current upper bound for the next interval
            lower_time_index = upper_time_index;
        }
    }
    return C_kj_matrix; // Return the computed regression weight matrix
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
double *SFS_theo(double *thetas, double *weight, int n, int m)
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
                n, 1, m, 1.0, weight, m, thetas, 1, 0.0, sfs_theo, 1);

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
double *SFS_to_freq(double *sfs_theo, int n)
{
    double n_snp = 0;                               // Variable to hold the total number of SNPs
    double *sfs_theo2 = malloc(sizeof(double) * n); // Allocate memory for the frequency array
    if (sfs_theo2 == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    // Calculate the total number of SNPs by summing the counts in sfs_theo
    for (int i = 0; i < n; i++)
        n_snp += sfs_theo[i];
    // Compute the frequency for each site and store in sfs_theo2
    for (int i = 0; i < n; i++)
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
double log_likelihood(double *sfs_obs, double *sfs_theo, int n)
{
    double llikelihood = 0;                       // Variable to hold the calculated log likelihood
    double *sfs_theo2 = SFS_to_freq(sfs_theo, n); // Convert theoretical SFS to frequencies
    // Calculate the log likelihood
    for (int i = 0; i < n; i++)
        llikelihood += (sfs_obs[i] * log(sfs_theo2[i]));

    free(sfs_theo2);    // Free the memory allocated for the frequency SFS
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
void replace_negative_with_1(double *theta, size_t n)
{
    for (size_t i = 0; i < n; i++) // Iterate through each element of the vector.
    {
        if (theta[i] < 0) // Check if the current value is negative.
        {
            theta[i] = 1.; // Replace the negative value with 1.
        }
    }
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
void system_resolution(solution *sol, double **sfs, double **cumul_weight, int sfs_length)
{
    // Step 1: Assemble regression weights from cumulative branch lengths
    double *weight = weight_assembly_1d(sfs_length, sol->nb_breakpoints + 1, sol->breakpoints, cumul_weight);
    double *regressors = calloc((sol->nb_breakpoints + 1) * sfs_length, sizeof(double));
    if (regressors == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    // Step 3: Compute the regression matrix using the regression weights (X^T X)-1X^T
    regressor_matrix(sfs_length, sol->nb_breakpoints + 1, weight, regressors);
    // Step 4: Estimate population mutation rates (thetas) using matrix multiplication (X^T X)-1X^T * SFS
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                sol->nb_breakpoints + 1, 1, sfs_length, 1.0, regressors, sfs_length, sfs[0], 1, 0.0, sol->thetas, 1);
    replace_negative_with_1(sol->thetas, sol->nb_breakpoints + 1); // if thetas are negatives they are replaced by 1 as population sizes cannot be inferior to 0
    double *sfs_theo = SFS_theo(sol->thetas, weight, sfs_length, sol->nb_breakpoints + 1);
    sol->log_likelihood = log_likelihood(sfs[1], sfs_theo, sfs_length);
    sol->distance = distance(sfs[0], sfs_theo, sfs_length);
    free(sfs_theo);
    free(weight);
    free(regressors);

    
}

/**
 * Generates the next combination of breakpoints in lexicographic order within a specified grid.
 *
 * @param brk            Array representing the current combination of breakpoints. breakpoints are index of times of changes from which cumulative branch lengthes have been computed
 *                       This array will be modified in place to hold the next combination.
 * @param nb_breakpoints The number of breakpoints (i.e., the size of the `brk` array).
 * @param grid_size      The size of the grid, representing the maximum index range for breakpoints.
 *
 * @return               1 if there are more combinations to generate,
 *                       0 if the last combination has been reached.
 */
int recursive_bk_combination(int *brk, int nb_breakpoints, int grid_size)
{
    // Start with the last breakpoint in the array
    int i = nb_breakpoints - 1;
    // Attempt to increment the last breakpoint by 1
    brk[i] += GRIDREFINE;
    // Adjust previous breakpoints if the current one exceeds allowable range
    while (i >= 0 && brk[i] / GRIDREFINE > grid_size + 1 - nb_breakpoints + i)
    {
        // Move to the previous breakpoint in the array
        i--;
        // Increment the previous breakpoint
        brk[i] += GRIDREFINE;
    }
    // For any indices after `i`, reset them to maintain increasing order
    for (int j = i + 1; j < nb_breakpoints; j++)
        brk[j] = brk[j - 1] + GRIDREFINE;
    // If the first breakpoint reaches its upper limit, all combinations have been generated
    if (i == 0 && brk[0] / GRIDREFINE == grid_size + 1 - nb_breakpoints)
        return 0;
    // Return 1 to indicate more combinations are available
    return 1;
}

// /**
//  * Converts relative times in the `time` array to absolute times.
//  *
//  * Relative times are initially derived from time indices on a grid and represent intervals scaled by
//  * effective population sizes at each historical breakpoint. This function adjusts each relative time to
//  * an absolute timescale, where times are expressed in units of Ne (effective population size) multiplied
//  * by the number of generations.
//  *
//  * Each interval is scaled by the population size relative to the present, converting times in `time`
//  * into cumulative (absolute) times. The absolute times represent the total elapsed time since the present.
//  *
//  * @param time           Array of relative times, modified in place to store absolute times.
//  * @param thetas         Array of population mutation rates, with each rate proportional to the Ne
//  *                        at its corresponding time interval.
//  * @param nb_breakpoints Number of breakpoints, representing the population size change intervals.
//  */
// void relatif_to_absolute_time(double *time, double *thetas, int nb_breakpoints)
// {
//     // Adjust relative times to absolute by scaling each by the relative Ne and accumulating past times.
//     for (int i = nb_breakpoints - 1; i >= 1; i--)
//         time[i] = (time[i] - time[i - 1]) * thetas[i] / thetas[0];
//     // Accumulate times to transform each into an absolute time from the present.
//     for (int i = 1; i < nb_breakpoints; i++)
//         time[i] += time[i - 1];
// }

/**
 * Checks if all mutation rate values (thetas) are positive.
 *
 * This validation ensures that all effective population size estimates, which are proportional to thetas,
 * are positive values. A negative or zero value would indicate an invalid population size.
 *
 * @param thetas         Array of population mutation rates.
 * @param nb_breakpoints Number of breakpoints (size of the `thetas` array minus one).
 *
 * @return               1 if all values are positive; 0 otherwise.
 */
int all_positive(double *thetas, int nb_breakpoints)
{
    for (int i = 0; i <= nb_breakpoints; i++)
    {
        if (thetas[i] <= 0.)
            return 0; // Found a non-positive value.
    }
    return 1; // All values are positive.
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
void save_solution(solution sol, int n_sample, char *out_file, double const_ren, int sfs_length)
{
    // Ensure the directory exists before proceeding
    FILE *file = fopen(out_file, "a");
    if (file == NULL)
    {
        perror("Error opening file");
        free(out_file);
        return;
    }

    // Write solution data to file
    fprintf(file, " >\n log_likelihood: %f ", sol.log_likelihood);
    fprintf(file, "\n distance: %f ", sol.distance);
    fprintf(file, "\n breakpoints: ");
    for (int i = 0; i < sol.nb_breakpoints; i++)
        fprintf(file, "%d ", sol.breakpoints[i]);
    fprintf(file, "\n thetas: ");
    for (int i = 0; i <= sol.nb_breakpoints; i++)
        fprintf(file, "%f ", sol.thetas[i] * const_ren);
    if (sol.se_thetas != NULL)
    {
        fprintf(file, "\n se_thetas: ");
        for (int i = 0; i <= sol.nb_breakpoints; i++)
            fprintf(file, "%f ", sol.se_thetas[i] * const_ren);
    }
    if (sol.residues != NULL)
    {
        fprintf(file, "\n residues: ");
        for (int i = 0; i < sfs_length; i++)
            fprintf(file, "%f ", sol.residues[i] * const_ren);
    }
    if (sol.fitted_sfs != NULL)
    {
        fprintf(file, "\n fitted_sfs: ");
        for (int i = 0; i < sfs_length; i++)
            fprintf(file, "%f ", sol.fitted_sfs[i] * const_ren);
    }
    fprintf(file, "\n");
    fclose(file);
    free(sol.residues);
    free(sol.se_thetas);
    free(sol.fitted_sfs);
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
void thetas_se(solution *sol, int sfs_length, double **cumul_weight)
{
    int nb_thetas = sol->nb_breakpoints + 1;

    // Calculate error variance using the residual sum of squares
    double error_var = sol->distance * sol->distance / (sfs_length - nb_thetas);

    // Assemble regression weights
    double *weight = weight_assembly_1d(sfs_length, nb_thetas, sol->breakpoints, cumul_weight);

    // Allocate memory for variance-covariance matrix
    double *XXT = calloc(nb_thetas * nb_thetas, sizeof(double));
    sol->se_thetas = calloc(nb_thetas, sizeof(double));

    // Compute variance-covariance matrix
    varcovar(sfs_length, nb_thetas, weight, XXT);

    // Calculate standard errors for each theta
    for (int j = 0; j < nb_thetas; j++)
        sol->se_thetas[j] = sqrt(error_var * XXT[j * nb_thetas + j]);

    // Free allocated memory
    free(weight);
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
void residues(solution *sol, double **cumul_weight, int sfs_length, double *sfs)
{
    int nb_thetas = sol->nb_breakpoints + 1;
    double *weight = weight_assembly_1d(sfs_length, nb_thetas, sol->breakpoints, cumul_weight);
    sol->fitted_sfs = SFS_theo(sol->thetas, weight, sfs_length, sol->nb_breakpoints + 1);
    sol->residues = calloc(sfs_length, sizeof(double));
    for (int i = 0; i < sfs_length; i++)
        sol->residues[i] = sfs[i] - sol->fitted_sfs[i];
    free(weight);
}



/**
 * Generates and evaluates all possible breakpoint combinations to find the optimal solution
 * with the highest log-likelihood, saving each intermediate solution and returning the best one.
 *
 * @param nb_breakpoints  Number of breakpoints (population size changes) to consider.
 * @param sfs_length      Length of the Site Frequency Spectrum (SFS).
 * @param cumul_weight    Cumulative weights matrix for branch lengths across time intervals.
 * @param sfs             Observed SFS data (theoretical and observed SFS).
 * @param grid_size       Total number of grid points available for placing breakpoints.
 * @param n_sample        Sample identifier or index used for saving solutions.
 *
 * @return                The `solution` structure containing the optimal combination of breakpoints
 *                        with the highest log-likelihood value.
 */
solution generate_brk_combinations(int nb_breakpoints, int sfs_length, double **cumul_weight, double **sfs, int grid_size, int n_sample)
{
    // Initialize the solution with specified number of breakpoints
    solution sol = init_solution_size(nb_breakpoints);
    // Set the last breakpoint to the grid limit (end boundary)
    sol.breakpoints[nb_breakpoints] = GRIDREFINE * grid_size + 1;
    // Flag for stopping the combination generation
    int arret = 1;
    // Calculate the initial solution with the given breakpoints
    system_resolution(&sol, sfs, cumul_weight, sfs_length);
    // Initialize a temporary solution to keep track of the best found solution
    solution tmp_sol = copy_solution(sol);
    // Loop through all breakpoint combinations until `arret` is set to 0
    while (arret && nb_breakpoints > 0)
    {
        // Generate the next combination of breakpoints
        arret = recursive_bk_combination(sol.breakpoints, nb_breakpoints, grid_size);
        // Calculate the solution (log-likelihood and distance) for the new combination
        system_resolution(&sol, sfs, cumul_weight, sfs_length);
        // If the new solution has a higher log-likelihood and all theta values are positive, keep its
        if (isnan(tmp_sol.log_likelihood) || sol.log_likelihood > tmp_sol.log_likelihood) // || (all_positive(sol.thetas, sol.nb_breakpoints) && !all_positive(tmp_sol.thetas, sol.nb_breakpoints))) //
        {
            //if (!all_positive(sol.thetas, sol.nb_breakpoints) && all_positive(tmp_sol.thetas, sol.nb_breakpoints))
                //continue;
            clear_solution(tmp_sol);
            tmp_sol = copy_solution(sol);
        }
    }
    clear_solution(sol);
    thetas_se(&tmp_sol, sfs_length, cumul_weight);
    residues(&tmp_sol, cumul_weight, sfs_length, sfs[0]);
    // Return the best solution with the optimal breakpoint combination for a givent nulber of changes in population size
    return tmp_sol;
}



int check_new(solution sol, solution sol_initiale, int breakpoint)
{
    if(sol.breakpoints[breakpoint] > sol.breakpoints[breakpoint + 1] || sol.breakpoints[breakpoint] < sol.breakpoints[breakpoint - 1])
        return 0;
    return 1;
}

solution refine_solution_b(solution sol_initiale, int b, double **sfs, double **cumul_weight, int sfs_length, int sign)
{
    solution sol1 = copy_solution(sol_initiale);
    solution solm = copy_solution(sol_initiale);
    sol1.breakpoints[b]  += sign;
    system_resolution(&sol1, sfs, cumul_weight, sfs_length);
    
    while(sol1.log_likelihood > solm.log_likelihood){
        solm = copy_solution(sol1);
        sol1.breakpoints[b]  += sign;
        system_resolution(&sol1, sfs, cumul_weight, sfs_length);
    }
    clear_solution(sol1);
    // clear_solution(sol2);
    return solm; 
}


void refine_solution(solution * sol_initiale, double **sfs, double **cumul_weight, int sfs_length)
{
   

    for(int i =0; i < 100; i ++)
    {   
        for(int b = sol_initiale->nb_breakpoints - 1; b >= 0; b --)
        {
            solution solp = refine_solution_b(*sol_initiale, b, sfs, cumul_weight, sfs_length, -1);
            solution solm = refine_solution_b(*sol_initiale, b, sfs, cumul_weight, sfs_length, +1);
            if(sol_initiale->log_likelihood < solm.log_likelihood)
            {
                clear_solution(*sol_initiale);
                *sol_initiale = copy_solution(solm);
            }
            else{
                if(sol_initiale-> log_likelihood < solp.log_likelihood)
                {
                    clear_solution(*sol_initiale);
                    *sol_initiale = copy_solution(solp);
                }
            }
            clear_solution(solp);
            clear_solution(solm);
        }
    }
    thetas_se(sol_initiale, sfs_length, cumul_weight);
    residues(sol_initiale, cumul_weight, sfs_length, sfs[0]);
}

/**
 * Finds the best scenarios by generating combinations of breakpoints and solving
 * the regression for each configuration to estimate population mutation rates (thetas).
 *
 * @param sfs_length   Length of the observed Site Frequency Spectrum (SFS).
 * @param cumul_weight Pointer to a 2D array of cumulative weights for regression.
 * @param sfs          Pointer to a 2D array where:
 *                      - sfs[0] is the training SFS used for regression.
 *                      - sfs[1] is the test SFS used for validation.
 * @param grid_size    size of the discrete gride of time points
 * @param n_sample     sample size
 * @param changes      Maximum number of allowed breakpoint changes.
 *
 * @return             Array of solutions containing the best solution for each number of 
 *                     population size changes, from 0 to `changes`. Each solution includes
 *                     estimated thetas and associated metrics (log-likelihood, distance).
 *
 * Steps:
 * 1. Calculate `const_ren`, the proportion of sites represented in the training SFS (`sfs[0]`).
 *    This is used to adjust scaling if needed.
 * 3. Iteratively generate scenarios for `nb_breakpoints` from 0 to `changes`:
 *    - Use `generate_brk_combinations` to compute the best solution for each number of breakpoints.
 * 4. Return the array of solutions for all configurations.
 */
solution *find_scenario(int sfs_length, double **cumul_weight, double **sfs, int grid_size, int n_sample, int changes)
{
    // Calculate the proportion of sites in the training SFS relative to the total
    double const_ren = (sfs[0][0] + sfs[1][0]) / sfs[0][0];
    if ((int)const_ren == 2) // If the value is exactly 2, set it to 1
        const_ren = 1.;

    // Allocate memory for storing solutions for all configurations of breakpoints
    solution *liste_solution = malloc(sizeof(solution) * (changes + 1));

    // Iterate through the number of breakpoints from 0 to `changes`
    int nb_breakpoints = 0;
    while (nb_breakpoints <= changes)
    {
        // Generate the best solution for the current number of breakpoints
        liste_solution[nb_breakpoints] = generate_brk_combinations(nb_breakpoints, sfs_length, cumul_weight, sfs, grid_size, n_sample);
        if (nb_breakpoints >= 1)
        {
            refine_solution(&liste_solution[nb_breakpoints], sfs, cumul_weight, sfs_length);
            // system_resolution(&liste_solution[nb_breakpoints], sfs, cumul_weight, sfs_length);
        }
        nb_breakpoints++;
    }

    // Return the list of solutions
    return liste_solution;
}


void copy_and_insert_in_sorted_array(int *breakpoints, int nb_breakpoint, int new_breakpoint, int *new_breakpoints)
{
    int i, j;
    int inserted = 0;

    // Parcours du tableau original et insertion dans le nouveau tableau
    for (i = 0, j = 0; i < nb_breakpoint; i++, j++)
    {
        // Insérer le nouvel élément à la bonne position
        if (!inserted && breakpoints[i] > new_breakpoint)
        {
            new_breakpoints[j] = new_breakpoint;
            inserted = 1;
            j++; // Pour éviter d'écraser l'élément suivant
        }
        // Copier l'élément du tableau original
        new_breakpoints[j] = breakpoints[i];
    }
}

solution * recent_infrence(solution *list_solution, int changes, double **sfs, double ** cumul_weight, int sfs_length, int n_sample)
{
    solution *list_solution_r = malloc(sizeof(solution) * (changes + 1));
    int br = 1;
    for(int i = 0; i <= changes; i++)
    {
        list_solution_r[i] = init_solution_size(list_solution[i].nb_breakpoints + 1);
        copy_and_insert_in_sorted_array(list_solution[i].breakpoints, list_solution_r[i].nb_breakpoints, br, list_solution_r[i].breakpoints);
        if(list_solution[i].breakpoints[0] != 1)    
            system_resolution(&list_solution_r[i], sfs, cumul_weight, sfs_length);
        else
        {
            list_solution_r[i] = list_solution[i];
        }
        thetas_se(&list_solution_r[i], sfs_length, cumul_weight);
        residues(&list_solution_r[i], cumul_weight, sfs_length, sfs[0]); 
        //save_solution(list_solution_r[i], n_sample, "a", out_file, const_ren, sfs_length);
    } 
    return list_solution_r;
}

// Fonction de comparaison pour qsort
int comparer(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

// Fonction pour compter les indices inférieurs à une valeur donnée
double inferior_count(int *tirages, int taille, double valeur)
{
    double count = 0;
    for (int i = 0; i < taille; i++)
    {
        if ((double)tirages[i] < valeur)
        {
            count += 1.;
        }
    }
    return count;
}

// Fonction pour tirer un tiers des indices sans remise et les trier
void site_sampling(int n, int *tirages, int taille)
{
    int *indices = malloc(sizeof(int) * n); // Tableau pour stocker les indices de 0 à n-1
    // Initialiser le tableau des indices de 0 à n-1
    for (int i = 0; i < n; i++)
        indices[i] = i;
    // Initialiser le générateur de nombres aléatoires
    srand(time(NULL));
    // Tirer au hasard un tiers des indices sans remise
    for (int i = 0; i < taille; i++)
    {
        int randIndex = rand() % (n - i);        // Tirer un indice aléatoire dans les indices restants
        tirages[i] = indices[randIndex];         // Stocker l'indice tiré
        indices[randIndex] = indices[n - i - 1]; // Remplacer l'indice tiré par le dernier non tiré
    }
    // Trier la liste des tirages
    qsort(tirages, taille, sizeof(int), comparer);
}

// Fonction pour calculer la somme cumulée des éléments de sfs
void cumulate_sfs(double *sfs, int taille, double *sommeCumulee)
{
    sommeCumulee[0] = sfs[0];
    for (int i = 1; i < taille; i++)
        sommeCumulee[i] = sommeCumulee[i - 1] + sfs[i];
}

// Fonction pour calculer la différence de comptage entre indices
void difference_count(int *tirages, int taille_tirages, double *sommeCumulee, int size, double *count)
{
    for (int i = 0; i < size; i++)
        count[i] = inferior_count(tirages, taille_tirages, sommeCumulee[i]);
}

// Fonction pour mettre à jour le tableau count
void mettre_a_jour_count(double *count, int taille)
{
    for (int i = taille - 1; i >= 1; i--)
    {
        count[i] = count[i] - count[i - 1];
    }
}

double *test_split(double *sfs, int size, int frac)
{
    double n_sites = 0;
    for (int i = 0; i < size; i++)
        n_sites += sfs[i];
    int taille_tirages = n_sites / frac;
    int *tirages = malloc(sizeof(int) * taille_tirages);
    double *cumulated_sfs = malloc(sizeof(double) * size);
    double *sfs_test = malloc(sizeof(double) * size);
    if (frac == 1)
    {
        for (int i = 0; i < size; i++)
            sfs_test[i] = sfs[i];
        return sfs_test;
    }
    site_sampling(n_sites, tirages, taille_tirages);
    cumulate_sfs(sfs, size, cumulated_sfs);
    difference_count(tirages, taille_tirages, cumulated_sfs, size, sfs_test);
    mettre_a_jour_count(sfs_test, size);
    for (int i = 0; i < size; i++)
        sfs[i] -= sfs_test[i];
    free(tirages);
    free(cumulated_sfs);
    return sfs_test;
}

/**
 * Modifies the system for analyses when the SFS is folded. This occurs when alleles are not oriented.
 * The folding process combines the first half of the SFS with the second half and adjusts the cumulative weights accordingly.
 *
 * @param sfs               Pointer to a 2D array where:
 *                           - sfs[0] is the observed Site Frequency Spectrum (SFS).
 * @param cumulative_weight Pointer to a 2D array of cumulative weights used in regression.
 * @param sfs_length        Length of the SFS, representing the number of observations.
 * @param grid_size         Size of the grid for cumulative weights.
 *
 * This function folds the SFS by adding the second half of the SFS to the first half and adjusting the cumulative weights
 * accordingly to account for the absence of allele orientation.
 */
void fold_sfs(double **sfs, double **cumulative_weight, int sfs_length, int grid_size)
{
    for (int i = 0; i < sfs_length / 2; i++)
    {
        // Combine the corresponding elements from the first and second halves of the SFS
        sfs[0][i] += sfs[0][sfs_length - 1 - i];

        // Adjust cumulative weights accordingly
        for (int j = 0; j < grid_size + 2; j++)
            cumulative_weight[i][j] += cumulative_weight[sfs_length - 1 - i][j];
    }
}


/**
 * Modifies the system for analyses when the SFS is folded. This occurs when alleles are not oriented.
 * The folding process combines the first half of the SFS with the second half and adjusts the cumulative weights accordingly.
 *
 * @param sfs               Pointer to a 2D array where:
 *                           - sfs[0] is the observed Site Frequency Spectrum (SFS).
 * @param cumulative_weight Pointer to a 2D array of cumulative weights used in regression.
 * @param sfs_length        Length of the SFS, representing the number of observations.
 * @param grid_size         Size of the grid for cumulative weights.
 *
 * This function folds the SFS by adding the second half of the SFS to the first half and adjusting the cumulative weights
 * accordingly to account for the absence of allele orientation.
 */
void sigleton_ignore(double **sfs, double **cumulative_weight, int grid_size)
{
    // Combine the corresponding elements from the first and second halves of the SFS
    sfs[0][0] = 0.;
    // Adjust cumulative weights accordingly
    for (int j = 0; j < grid_size + 2; j++)
        cumulative_weight[0][j] = 0.;
}

/**
 * Modifies the system for analyses when the SFS is folded. This occurs when alleles are not oriented.
 * The folding process combines the first half of the SFS with the second half and adjusts the cumulative weights accordingly.
 *
 * @param sfs               Pointer to a 2D array where:
 *                           - sfs[0] is the observed Site Frequency Spectrum (SFS).
 * @param cumulative_weight Pointer to a 2D array of cumulative weights used in regression.
 * @param sfs_length        Length of the SFS, representing the number of observations.
 * @param grid_size         Size of the grid for cumulative weights.
 *
 * This function folds the SFS by adding the second half of the SFS to the first half and adjusting the cumulative weights
 * accordingly to account for the absence of allele orientation.
 */
void singleton_erased(double **sfs, double **cumulative_weight, int sfs_length)
{
    for (int i = 0; i < sfs_length - 1; i++)
    {
        // Combine the corresponding elements from the first and second halves of the SFS
        sfs[0][i] = sfs[0][i + 1];
        // Adjust cumulative weights accordingly
        cumulative_weight[i] = cumulative_weight[i + 1];
    }
}

