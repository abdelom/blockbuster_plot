#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "blockbuster_grid.h"


// void free_integral_grid(double **cjk_grid, int n_sample)
// {
//     for (int i = 0; i < n_sample - 1; i++)
//         free(cjk_grid[i]); // Free each sub-array
//     free(cjk_grid);        // Free the main array
// }

void matrix_multiply(int n, int t, double *A, double **B, double **C)
{
    /**
     * @brief Performs matrix multiplication and stores the result in a given matrix.
     *
     * This function multiplies a vector A (1D array) by a matrix B (2D array) and
     * stores the result in matrix C (2D array). The computation is done for a specified
     * number of rows and columns, and the final results in matrix C are adjusted by
     * subtracting the inverse of the row index and rounding to six decimal places.
     *
     * @param n The number of rows in the output matrix C, which corresponds to n_sample - 1.
     * @param t The number of columns in matrix B, indicating the number of output columns.
     * @param A A pointer to a 1D array representing the first matrix (vector).
     * @param B A pointer to a 2D array representing the second matrix.
     * @param C A pointer to a 2D array where the result of the multiplication will be stored.
     */
    printf("%d \n", t); // Print the number of columns in matrix B
    for (int j = 0; j < t - 2; j++)
    { // Loop over each column of matrix B (except the last two)
        for (int i = 0; i < (n - 1); i++)
        {                                                              // Loop over each row of output matrix C
            for (int k = 0; k < (n - 1); k++)                          // Loop to perform the multiplication
                C[i][j + 1] += A[i * (n - 1) + k] * B[(n - 2) - k][j]; // Multiply A and B, accumulate in C
            // if (C[i][j] < 1e-6){ // Commented out: conditional adjustment
            //     C[i][j] = 0; // Set small values to zero
            //     break; // Break the loop if a condition is met
            C[i][j + 1] -= 1. / (double)(i + 1);           // Adjust C by subtracting the inverse of the index
            C[i][j + 1] = -C[i][j + 1];                    // Negate the value in C
            C[i][j + 1] = round(C[i][j + 1] * 1e6) * 1e-6; // Round to six decimal places
        }
    }
}


long double *init_wik(int n_sample, int sfs_length)
{
    long double *W_i_k = calloc((n_sample - 1) * sfs_length, sizeof(long double)); // Allocate memory for the array
    if (!W_i_k)                                                                        // Check for successful allocation
        exit(3);
    for (int i = 1; i <= sfs_length; i++)
    {
        W_i_k[(i - 1) * (n_sample - 1)] = 6. / (long double)(n_sample + 1);
        W_i_k[(i - 1) * (n_sample - 1) + 1] = 30. * (long double)(n_sample - 2 * i) / (long double)((n_sample + 1) * (n_sample + 2));
    }
    return W_i_k;
}


/**
 * Generates a list of numbers logarithmically spaced between two given bounds, writes them to a file,
 * and returns them in a dynamically allocated array.
 *
 * @param min      The lower bound of the logarithmic scale (must be > 0).
 * @param max      The upper bound of the logarithmic scale (must be > min).
 * @param n        The number of points to generate (must be > 0).
 * @return         A pointer to the dynamically allocated array of points, or NULL if an error occurs.
 */
// Function to generate a logarithmic scale with n+1 points, ensuring the last point is 2 if 2 < max.
double* generate_logarithmic_scale(int grid_size, double upper_bound, double lower_bound) 
{
    // Calculate the actual number of points
    // Allocate memory for the points
    double* log_scale = malloc(grid_size * sizeof(double));
    if (!log_scale) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }    

    // Generate logarithmic points
    double log_min = log10(lower_bound);
    double log_max = log10(upper_bound);
    double step = (log_max - log_min) / (grid_size - 1);

    for (int i = 0; i < grid_size; i++)
        log_scale[i] = pow(10, log_min + i * step);

    
    // log_scale[grid_size - 1] = log_scale[grid_size - 2] + 1.0;
    return log_scale;
}


/**
 * Generates a list of numbers linearly spaced between two given bounds, writes them to a file,
 * and returns them in a dynamically allocated array.
 *
 * @param grid_size  The number of points to generate (must be > 0).
 * @param upper_bound The upper bound of the linear scale.
 * @param lower_bound The lower bound of the linear scale.
 * @param file_name  The name of the file where the points will be saved.
 * @return           A pointer to the dynamically allocated array of points, or NULL if an error occurs.
 */
double* generate_linear_scale(int grid_size, double upper_bound, double lower_bound) 
{
    // Allocate memory for the points
    double* linear_scale = malloc(grid_size * sizeof(double));
    if (!linear_scale) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    // Generate linearly spaced points
    double step = (upper_bound - lower_bound) / (grid_size - 2);

    for (int i = 0; i < grid_size - 1; i++) {
        linear_scale[i] = lower_bound + i * step;
    }

    // Ensure the last point is 2.0 (if it fits within the range)
    linear_scale[grid_size - 1] = 2.0;
    return linear_scale;
}

void save_cumulated_weight(int sfs_length, int grid_size, double **matrix, char *filename)
{
    /**
     * @brief Saves the cumulative weight matrix to a specified file.
     * This function writes a 2D array (matrix) representing cumulative branch lengths
     * to a file in a space-separated format. Each element of the matrix is written
     * with six decimal places. The function checks for errors during file opening
     * and prints a success message upon completion.
     * @param n_sample The sample size (number of individuals) corresponding to the number of rows in the matrix.
     * @param grid_size The size of the grid corresponding to the number of columns in the matrix.
     * @param matrix A pointer to a 2D array (double**) containing the cumulative weights to be saved.
     * @param filename A pointer to a string representing the name of the file where the matrix will be saved.
     */
    FILE *f = fopen(filename, "w"); // Open the specified file for writing
    if (f == NULL)
    {                                                           // Check if the file opened successfully
        printf("Erreur d'ouverture du fichier %s\n", filename); // Print an error message if not
        exit(1);                                                // Exit the program
    }
    for (int i = 0; i < sfs_length; i++)
    { // Loop over each row
        for (int j = 0; j < grid_size; j++)
        {                                      // Loop over each column
            fprintf(f, "%.6f ", matrix[i][j]); // Write each element with six decimal places
        }
        fprintf(f, "\n"); // New line after each row
    }

    fclose(f);                                                             // Close the file
    printf("Matrice écrite dans le fichier %s avec succès !\n", filename); // Print success message
}


long double *Wik(int n_sample, int sfs_length)
{
    long double Wk, Wkp1;
    long double tmp_k, tmp_kp1;
    long double *W_i_k = init_wik(n_sample, sfs_length);
    for (int i = 1; i <= sfs_length; i++)
    {
        // printf("%d\n", i);
        for (int k = 2; k <= n_sample - 2; k++)
        {
            Wk = W_i_k[(i - 1) * (n_sample - 1) + (k - 2)];
            Wkp1 = W_i_k[(i - 1) * (n_sample - 1) + (k - 1)];
            tmp_k = (long double)((1 + k) * (3 + 2 * k) * (n_sample - k)) / (long double)(k * (2 * k - 1) * (n_sample + k + 1));
            tmp_kp1 = (long double)((3 + 2 * k) * (n_sample - 2 * i)) / (long double)(k * (n_sample + k + 1));
            W_i_k[(i - 1) * (n_sample - 1) + k] = -tmp_k * Wk + tmp_kp1 * Wkp1;
        //     if (W_i_k[(i - 1) * (n_sample - 1) + k] < 1e-7 && W_i_k[(i - 1) * (n_sample - 1) + k] > -1e-7)
        //     {
        //         if (k % 2 == 0)
        //             break;
        //     }
        }
    }
    return W_i_k;
}


long double *element_time(int n_sample, double Hj)
{
    long double *vkj = malloc(sizeof(long double) * (n_sample - 1));
    double two_k;
    for (int k = 2; k <= n_sample; k++)
    {
        two_k = (double)(k * (k - 1) / 2);
        vkj[k - 2] = (1. - exp(-two_k * Hj)) / (2. * two_k);
    }
    return vkj;
}


void weigth_grid_i(double **weight_grid, long double *wik, int col, int n_sample, int sfs_length, double Hj)
{
    long double *vkj = element_time(n_sample, Hj);
    for (int j = 0; j < sfs_length; j++)
    {
        for (int k = 0; k < n_sample - 1; k++)
            weight_grid[j][col] += vkj[k] * wik[j * (n_sample - 1) + k];
            // if(wik[j * (n_sample - 1) + k] < 1e-8)
            // break;
    }
    free(vkj);
}


/**
 * Computes the cumulative branch lengths for a given set of sample sizes and time intervals
 * based on either a logarithmic or linear scale for time.
 * This function generates weights used for regression by computing the cumulative branch lengths
 * for different descendant groups within the specified time intervals.
 *
 * @param n_sample     The number of samples.
 * @param grid_size    The size of the grid representing the time intervals.
 * @param upper_bound  The upper bound for the time scale.
 * @param lower_bound  The lower bound for the time scale.
 * @param file_name    The file name for saving the scale (if needed).
 * @param log          Flag to choose between logarithmic (1) or linear (0) scale for time.
 *
 * @return A 2D array `weight_grid` of size `(n_sample - 1) * (grid_size + 2)` containing 
 *         the cumulative branch lengths used in regression, where each entry 
 *         `weight_grid[i][j]` represents the cumulative branch length for `i + 1` descendants
 *         in the present up to time interval `j`. The caller is responsible for freeing this memory.
 *
 * Example usage:
 *     double **weights = cumulatve_weight_v2(n_sample, grid_size, upper_bound, lower_bound, file_name, log);
 *     // Use `weights` here
 *     free(weights);  // Free the allocated memory when done
 */
double **cumulatve_weight_v2(int n_sample, int sfs_length, int grid_size, double *H)
{


    // See kimmel and polanski 2003
    long double *wik = Wik(n_sample, sfs_length); 

    // Allocate memory for the weight grid, which stores cumulative branch lengths
    double **weight_grid = malloc(sfs_length * sizeof(double *)); // (n_sample - 1) rows for different descendant groups
    for (int i = 0; i < sfs_length; i++)
        weight_grid[i] = calloc(grid_size + 2, sizeof(double)); // Initialize each row with size grid_size + 2

    // Calculate the cumulative branch lengths for each descendant group in the grid
    for (int i = 1; i < grid_size + 1; i++)
        weigth_grid_i(weight_grid, wik, i, n_sample, sfs_length, H[i - 1]);

    // Finalize the cumulative branch lengths for the last time interval (infinity)
    weigth_grid_i(weight_grid, wik, grid_size + 1, n_sample, sfs_length, INFINITY);

    // Free allocated memory for the time scale and P_i,k values
    // free(H);
    free(wik);

    // Return the computed cumulative weight grid
    return weight_grid;
}


void fold_time_grid(Time_gride *tg, SFS sfs)
{
    if(!sfs.oriented)
    { 
        int folded_size = (sfs.n_haplotypes - 1) / 2 + (sfs.n_haplotypes - 1) % 2;
        int unfolded_size = sfs.n_haplotypes - 1;
        if(sfs.troncation > 5)
        unfolded_size = (unfolded_size < sfs.troncation) ? unfolded_size : sfs.troncation; 
        for(int i = folded_size; i < unfolded_size; i ++)
        {
            for (int j = 0; j < tg->grid_size * 1000 + 3; j++)
            {
                if(!sfs.singleton && i == sfs.n_haplotypes - 2)
                    tg->cumulative_bl[sfs.n_haplotypes - 2 - i][j] = tg->cumulative_bl[i][j];
                else
                    tg->cumulative_bl[sfs.n_haplotypes - 2 - i][j] += tg->cumulative_bl[i][j];
            }
            free(tg->cumulative_bl[i]);
        }
    }
}


void erased_singleton_grid(Time_gride *tg, SFS sfs)
{
    if(sfs.oriented && !sfs.singleton)
    {
        double *tmp = tg->cumulative_bl[0];
        for(int i = 1; i < sfs.sfs_length; i ++)
        {
            tg->cumulative_bl[i - 1] = tg->cumulative_bl[i];
        }   
        free(tmp);
    }
}


Time_gride init_time_grid(SFS sfs, int grid_size, double ub, double lb)
{
    Time_gride time_grid;
    time_grid.grid_size = grid_size;
    time_grid.time_scale = generate_logarithmic_scale(grid_size * 1000 + 1, ub, lb);
    int sfs_length = sfs.n_haplotypes - 1;
    if(sfs.troncation > 5)
        sfs_length = (sfs_length < sfs.troncation) ? sfs_length : sfs.troncation;
    printf("%d\n", sfs_length);
    time_grid.cumulative_bl = cumulatve_weight_v2(
        sfs.n_haplotypes,
        sfs_length,
        grid_size * 1000 + 1,
        time_grid.time_scale
    );
    fold_time_grid(&time_grid, sfs);
    erased_singleton_grid(&time_grid, sfs);
    return time_grid;
}


Time_gride init_time_grid_H(int n_haplotypes, int grid_size, double *H)
{
    Time_gride time_grid;
    time_grid.grid_size = grid_size;
    time_grid.time_scale = H;
    time_grid.cumulative_bl =  cumulatve_weight_v2(n_haplotypes, n_haplotypes - 1, grid_size, H);
    return time_grid;
}


void clear_time_grid(Time_gride tg, int sfs_length)
{
    free(tg.time_scale);
    for (int i = 0; i < sfs_length; i++)
        free(tg.cumulative_bl[i]); // Free each sub-array
    free(tg.cumulative_bl);
}