#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "blockbuster_grid.h"

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


/**
 * @brief Initializes the W_i_k matrix used for cumulative branch length calculations.
 *
 * This function allocates and initializes a matrix of size (sfs_length x (n_sample - 1))
 * that contains constants required for computing the cumulative branch lengths
 * in the `Time_gride` structure.
 *
 * The first two entries of each row are initialized with specific formulas,
 * and the remaining entries are expected to be computed recursively later in `Wik()`.
 *
 * @param n_sample    Number of haplotypes (sample size).
 * @param sfs_length  Number of SFS bins.
 *
 * @return Pointer to the dynamically allocated and initialized W_i_k array.
 *         Exits the program with code 3 if memory allocation fails.
 */
double *init_wik(int n_sample, int sfs_length)
{
    // Allocate memory for the W_i_k array (rows = sfs_length, cols = n_sample - 1)
    double *W_i_k = calloc((n_sample - 1) * sfs_length, sizeof(double));
    if (!W_i_k)  // Check allocation success
        exit(3);

    // Initialize the first two entries of each row
    for (int i = 1; i <= sfs_length; i++)
    {
        W_i_k[(i - 1) * (n_sample - 1)] = 6.0L / (double)(n_sample + 1);
        W_i_k[(i - 1) * (n_sample - 1) + 1] =
            30.0L * (double)(n_sample - 2 * i) /
            (double)((n_sample + 1) * (n_sample + 2));
    }

    return W_i_k;
}

/**
 * @brief Generates a list of logarithmically spaced time points between two given bounds.
 *
 * This function creates a logarithmic time scale between `lower_bound` and `upper_bound`
 * and returns the resulting points in a dynamically allocated array.
 * The cumulative branch length matrix from the `Time_gride` structure
 * is computed on this logarithmic time scale.
 *
 * @param grid_size     Number of time points to generate (must be > 0).
 * @param upper_bound   The upper bound of the logarithmic scale (must be > lower_bound).
 * @param lower_bound   The lower bound of the logarithmic scale (must be > 0).
 *
 * @return A pointer to the dynamically allocated array containing the time points,
 *         or NULL if a memory allocation error occurs.
 *
 * @note The time scale is logarithmically spaced using base 10.
 *       The generated points can be used to compute cumulative branch lengths
 *       or other time-dependent quantities in coalescent-based models.
 */
double* generate_logarithmic_scale(int grid_size, double upper_bound, double lower_bound) 
{
    // Allocate memory for the logarithmic scale
    double* log_scale = malloc(grid_size * sizeof(double));
    if (!log_scale) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    // Compute logarithmic spacing parameters
    double log_min = log10(lower_bound);
    double log_max = log10(upper_bound);
    double step = (log_max - log_min) / (grid_size - 2);

    // Generate logarithmically spaced time points
    for (int i = 0; i < grid_size; i++)
        log_scale[i] = pow(10, log_min + i * step);

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


/**
 * @brief Calcule la matrice des constantes W(i, k) utilisées pour le calcul
 *        des longueurs de branches cumulées dans la structure `Time_gride`.
 *
 * Cette fonction génère une matrice de taille (sample_size - 1) × (sample_size - 1)
 * contenant les constantes nécessaires au calcul des longueurs de branches cumulées
 * dans le cadre de modèles de coalescence.  
 *
 * Chaque valeur W(i, k) est obtenue par récurrence à partir des valeurs précédentes,
 * selon des coefficients dépendant de la taille de l’échantillon et du nombre de descendants. voir kimmel et polanski 2003
 *
 * @param n_sample      Taille de l’échantillon (nombre d’individus simulés).
 * @param sfs_length    Longueur du SFS (nombre de classes de fréquence).
 *
 * @return Un pointeur vers la matrice W(i, k) allouée dynamiquement
 *         (de taille (sfs_length) × (n_sample - 1))
 *
 */
double *Wik(int n_sample, int sfs_length)
{
    double Wk, Wkp1;
    double tmp_k, tmp_kp1;
    double *W_i_k = init_wik(n_sample, sfs_length);

    for (int i = 1; i <= sfs_length; i++)
    {
        for (int k = 2; k <= n_sample - 2; k++)
        {
            Wk = W_i_k[(i - 1) * (n_sample - 1) + (k - 2)];
            Wkp1 = W_i_k[(i - 1) * (n_sample - 1) + (k - 1)];

            tmp_k = (double)((1 + k) * (3 + 2 * k) * (n_sample - k))
                   / (double)(k * (2 * k - 1) * (n_sample + k + 1));

            tmp_kp1 = (double)((3 + 2 * k) * (n_sample - 2 * i))
                     / (double)(k * (n_sample + k + 1));

            W_i_k[(i - 1) * (n_sample - 1) + k] = -tmp_k * Wk + tmp_kp1 * Wkp1;
        }
    }
    return W_i_k;
}


/**
 * @brief Computes the elementary branch lengths for a given time interval.
 *
 * This function calculates, for each number of lineages k, an integral term betweeb 0 and a time Hj under the neutral coalescent model.
 * The resulting vector is used to weight the W(i, k) constants in order to obtain
 * the cumulative branch lengths over the time interval.
 *
 * @param n_sample  Sample size (number of lineages at present).
 * @param Hj        Time (upper bound of the interval considered).
 *
 * @return A pointer to a dynamically allocated array of size (n_sample - 1)
 *         containing the elementary branch lengths for each k.
 *         The array must be freed by the user using `free()`.
 */
double *element_time(int n_sample, double Hj)
{
    double *vkj = malloc(sizeof(double) * (n_sample - 1));
    if (!vkj)
    {
        fprintf(stderr, "Memory allocation failed in element_time().\n");
        return NULL;
    }

    double two_k;
    for (int k = 2; k <= n_sample; k++)
    {
        two_k = (double)(k * (k - 1) / 2);
        vkj[k - 2] = (1.0 - exp(-two_k * Hj)) / (2.0 * two_k);
    }

    return vkj;
}


/**
 * @brief Computes cumulative branch lengths for a given column of the `Time_gride.cumulative_bl` matrix.
 *
 * This function combines the W(i, k) constants with the integral term
 * computed by `element_time()` to obtain the weighted sum corresponding to the cumulative
 * branch length between 0 and Hj (the j-th point of the logarithmic time grid).
 *
 * @param weight_grid  Matrix of cumulative branch lengths to fill.
 * @param wik          Array of W(i, k) constants.
 * @param col          Column index in the matrix to compute.
 * @param n_sample     Sample size (number of lineages at present).
 * @param sfs_length   Length of the SFS (number of frequency classes).
 * @param Hj           Time corresponding to the column (upper bound of interval).
 */
void weigth_grid_i(double **weight_grid, double *wik, int col, int n_sample, int sfs_length, double Hj)
{
    double *vkj = element_time(n_sample, Hj);
    if (!vkj)
        return;

    for (int j = 0; j < sfs_length; j++)
    {
        weight_grid[j][col] = 0.0;
        for (int k = 0; k < n_sample - 1; k++)
            weight_grid[j][col] += vkj[k] * wik[j * (n_sample - 1) + k];
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
void cumulatve_weight_v2(int n_sample, int sfs_length, int grid_size, Time_gride *tg)
{
    // Allouer la grille des poids (weight_grid) seulement si elle n'est pas déjà allouée
    if (tg->cumulative_bl == NULL) {
        tg->cumulative_bl  = malloc(sfs_length * sizeof(double *)); // (n_sample - 1) lignes pour différents groupes de descendants
        for (int i = 0; i < sfs_length; i++)
            tg->cumulative_bl [i] = calloc(grid_size + 2, sizeof(double));
    }
    // Calculate the cumulative branch lengths for each descendant group in the grid
    for (int i = 1; i < grid_size + 1; i++)
        weigth_grid_i(tg->cumulative_bl ,tg->wik, i, n_sample, sfs_length, tg->time_scale[i - 1]);

    // Finalize the cumulative branch lengths for the last time interval (infinity)
    weigth_grid_i(tg->cumulative_bl , tg->wik, grid_size + 1, n_sample, sfs_length, INFINITY);
}


/**
 * @brief Folds the cumulative branch length matrix according to the SFS orientation.
 *
 * If the SFS is unoriented, this function combines the rows of the cumulative branch
 * length matrix in a mirrored fashion: the first half is summed with the second half.
 * The last row is treated differently if singletons are ignored. Folded rows beyond the
 * effective size are freed to save memory.
 *
 * @param tg  Pointer to the Time_gride structure containing the cumulative branch length matrix.
 * @param sfs SFS structure providing orientation, truncation, and singleton information.
 */
void fold_time_grid(Time_gride *tg, SFS sfs)
{
    if (!sfs.oriented)
    { 
        int folded_size = (sfs.n_haplotypes - 1) / 2 + (sfs.n_haplotypes - 1) % 2;
        int unfolded_size = sfs.n_haplotypes - 1;

        if (sfs.troncation > 5)
            unfolded_size = (unfolded_size < sfs.troncation) ? unfolded_size : sfs.troncation; 

        for (int i = folded_size; i < unfolded_size; i++)
        {
            for (int j = 0; j < tg->grid_size * 1000 + 3; j++)
            {
                // Special handling for ignored singleton row
                if (!sfs.singleton && i == sfs.n_haplotypes - 2)
                    tg->cumulative_bl[sfs.n_haplotypes - 2 - i][j] = tg->cumulative_bl[i][j];
                else
                    tg->cumulative_bl[sfs.n_haplotypes - 2 - i][j] += tg->cumulative_bl[i][j];
            }
            free(tg->cumulative_bl[i]);  // Free mirrored row after folding
        }
    }
}

/**
 * @brief Removes the singleton row from the cumulative branch length matrix if ignored.
 *
 * If the SFS is oriented and singletons are ignored, this function shifts all rows
 * of the cumulative branch length matrix up by one, effectively removing the first row,
 * and frees the original first row.
 *
 * @param tg  Pointer to the Time_gride structure containing the cumulative branch length matrix.
 * @param sfs SFS structure providing orientation and singleton information.
 */
void erased_singleton_grid(Time_gride *tg, SFS sfs)
{
    if (sfs.oriented && !sfs.singleton)
    {
        double *tmp = tg->cumulative_bl[0];
        for (int i = 1; i < sfs.sfs_length; i++)
        {
            tg->cumulative_bl[i - 1] = tg->cumulative_bl[i];
        }
        free(tmp);  // Free the original singleton row
    }
}

/**
 * @brief Initializes a Time_gride structure based on the SFS and a logarithmic time scale.
 *
 * This function performs the following steps:
 * 1. Generates a logarithmic time scale between the given lower and upper bounds.
 * 2. Computes the cumulative branch length matrix (`cumulative_bl`) for each SFS bin 
 *    and sample size using `cumulatve_weight_v2`.
 * 3. Applies SFS-based preprocessing to the cumulative branch length matrix:
 *    - **Folding:** If the SFS is unoriented, the first half of the rows is combined
 *      with the mirrored second half of the matrix.
 *    - **Truncation:** Rows beyond the truncation bin are ignored or set to zero.
 *    - **Singleton removal:** The first row (corresponding to singletons) is removed
 *      or set to zero if singletons are ignored.
 *
 * Note: The folding, truncation, and singleton removal operations are applied **to the rows
 * of the cumulative branch length matrix**, not directly to the SFS vector.
 *
 * @param sfs        SFS structure containing sample information and preprocessing options.
 * @param grid_size  Number of points in the base time grid.
 * @param ub         Upper bound of the time scale.
 * @param lb         Lower bound of the time scale.
 *
 * @return A Time_gride structure with initialized time scale and preprocessed cumulative branch lengths.
 */
Time_gride init_time_grid(SFS sfs, int grid_size, double ub, double lb)
{
    Time_gride time_grid;
    time_grid.grid_size = grid_size;
    time_grid.cumulative_bl = NULL;
    // Generate logarithmic time scale with grid refinement
    time_grid.time_scale = generate_logarithmic_scale(grid_size * GRIDREFINE + 1, ub, lb);

    // Determine effective SFS length for cumulative branch length calculation
    int sfs_length = sfs.n_haplotypes - 1;
    if (sfs.troncation > 5)
        sfs_length = (sfs_length < sfs.troncation) ? sfs_length : sfs.troncation;

    printf("%d\n", sfs_length);
        // See kimmel and polanski 2003
    time_grid.wik = Wik(sfs.n_haplotypes, sfs_length); 
    // Compute cumulative branch length matrix for all SFS bins and time points
    cumulatve_weight_v2(
        sfs.n_haplotypes,
        sfs_length,
        grid_size * GRIDREFINE + 1,
        &time_grid
    );

    // Apply folding and singleton removal to the rows of cumulative_bl
    fold_time_grid(&time_grid, sfs);
    erased_singleton_grid(&time_grid, sfs);

    return time_grid;
}


/**
 * @brief Initializes a Time_gride structure using a given time vector H.
 *
 * This function creates a Time_gride object, sets its time scale to the provided
 * array H, and computes the cumulative branch length matrix for each SFS class
 * using `cumulatve_weight_v2()`.
 *
 * @param n_haplotypes  Number of haplotypes (sample size).
 * @param grid_size     Number of points in the time grid.
 * @param H             Array of time points.
 *
 * @return A Time_gride structure with initialized time scale and cumulative branch lengths.
 */
Time_gride init_time_grid_H(int n_haplotypes, int grid_size, double *H)
{
    Time_gride time_grid;
    time_grid.grid_size = grid_size;
    time_grid.cumulative_bl = NULL;
    time_grid.time_scale = H;
    time_grid.wik = Wik(n_haplotypes, n_haplotypes - 1);
    cumulatve_weight_v2(n_haplotypes, n_haplotypes - 1, grid_size, &time_grid);
    return time_grid;
}


/**
 * @brief Initializes a Time_gride structure using a given time vector H.
 *
 * This function creates a Time_gride object, sets its time scale to the provided
 * array H, and computes the cumulative branch length matrix for each SFS class
 * using `cumulatve_weight_v2()`.
 *
 * @param n_haplotypes  Number of haplotypes (sample size).
 * @param grid_size     Number of points in the time grid.
 * @param H             Array of time points.
 *
 * @return A Time_gride structure with initialized time scale and cumulative branch lengths.
 */
Time_gride init_time_grid_H_wik(Time_gride tg, SFS sfs,  double *H)
{
    Time_gride time_grid;
    time_grid.grid_size = tg.grid_size;
    time_grid.cumulative_bl = NULL;
    time_grid.time_scale = H;
    time_grid.wik = tg.wik;
    cumulatve_weight_v2(sfs.n_haplotypes, sfs.n_haplotypes - 1, 6, &time_grid);
    clear_time_grid(tg, sfs.n_haplotypes - 1);
        // Apply folding and singleton removal to the rows of cumulative_bl
    fold_time_grid(&time_grid, sfs);
    erased_singleton_grid(&time_grid, sfs);
    return time_grid;
}

/**
 * @brief Frees the memory associated with a Time_gride structure.
 *
 * This function frees the dynamically allocated time scale array and each row
 * of the cumulative branch length matrix, followed by the matrix itself.
 *
 * @param tg           The Time_gride object to clear.
 * @param sfs_length   Number of SFS classes (number of rows in cumulative_bl).
 */
void clear_time_grid(Time_gride tg, int sfs_length)
{
    free(tg.time_scale);

    for (int i = 0; i < sfs_length; i++)
        free(tg.cumulative_bl[i]);  // Free each row of the cumulative branch length matrix

    free(tg.cumulative_bl);         // Free the matrix itself
}
