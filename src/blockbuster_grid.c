#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

void free_integral_grid(double **cjk_grid, int n_sample)
{
    for (int i = 0; i < n_sample - 1; i++)
        free(cjk_grid[i]); // Free each sub-array
    free(cjk_grid);        // Free the main array
}

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

long double *init_wik(int n_sample)
{
    long double *W_i_k = calloc((n_sample - 1) * (n_sample - 1), sizeof(long double)); // Allocate memory for the array
    if (!W_i_k)                                                                        // Check for successful allocation
        exit(3);
    for (int i = 1; i < n_sample; i++)
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
 * @param filename The name of the file where the points will be saved.
 * @param min      The lower bound of the logarithmic scale (must be > 0).
 * @param max      The upper bound of the logarithmic scale (must be > min).
 * @param n        The number of points to generate (must be > 0).
 * @return         A pointer to the dynamically allocated array of points, or NULL if an error occurs.
 */
// Function to generate a logarithmic scale with n+1 points, ensuring the last point is 2 if 2 < max.
double* generate_logarithmic_scale(int grid_size, double upper_bound, double lower_bound, char *file_name) 
{
    // Calculate the actual number of points
    // Allocate memory for the points
    double* log_scale = malloc(grid_size * sizeof(double));
    if (!log_scale) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }    
    // Open the file for writing
    FILE* file = fopen(file_name, "w");
    if (!file) {
        fprintf(stderr, "Failed to open the file.\n");
        free(log_scale);
        return NULL;
    }

    // Generate logarithmic points
    double log_min = log10(lower_bound);
    double log_max = log10(upper_bound);
    double step = (log_max - log_min) / (grid_size - 2);

    for (int i = 0; i < grid_size - 1; i++) {
        log_scale[i] = pow(10, log_min + i * step);
        fprintf(file, "%.10f\n", log_scale[i]);
    }

    
    log_scale[grid_size - 1] = 2.0;
    fprintf(file, "%.10f\n", log_scale[grid_size - 1]);
    fclose(file);
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
double* generate_linear_scale(int grid_size, double upper_bound, double lower_bound, char *file_name) 
{
    // Allocate memory for the points
    double* linear_scale = malloc(grid_size * sizeof(double));
    if (!linear_scale) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    // Open the file for writing
    FILE* file = fopen(file_name, "w");
    if (!file) {
        fprintf(stderr, "Failed to open the file.\n");
        free(linear_scale);
        return NULL;
    }

    // Generate linearly spaced points
    double step = (upper_bound - lower_bound) / (grid_size - 2);

    for (int i = 0; i < grid_size - 1; i++) {
        linear_scale[i] = lower_bound + i * step;
        fprintf(file, "%.10f\n", linear_scale[i]);
    }

    // Ensure the last point is 2.0 (if it fits within the range)
    linear_scale[grid_size - 1] = 2.0;
    fprintf(file, "%.10f\n", linear_scale[grid_size - 1]);

    fclose(file);
    return linear_scale;
}

void save_cumulated_weight(int n_sample, int grid_size, double **matrix, char *filename)
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

    for (int i = 0; i < n_sample - 1; i++)
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

long double *Wik(int n_sample)
{
    long double Wk, Wkp1;
    long double tmp_k, tmp_kp1;
    long double *W_i_k = init_wik(n_sample);
    for (int i = 1; i < n_sample; i++)
    {
        // printf("%d\n", i);
        for (int k = 2; k <= n_sample - 2; k++)
        {
            Wk = W_i_k[(i - 1) * (n_sample - 1) + (k - 2)];
            Wkp1 = W_i_k[(i - 1) * (n_sample - 1) + (k - 1)];
            tmp_k = (long double)((1 + k) * (3 + 2 * k) * (n_sample - k)) / (long double)(k * (2 * k - 1) * (n_sample + k + 1));
            tmp_kp1 = (long double)((3 + 2 * k) * (n_sample - 2 * i)) / (long double)(k * (n_sample + k + 1));
            W_i_k[(i - 1) * (n_sample - 1) + k] = -tmp_k * Wk + tmp_kp1 * Wkp1;
            if (W_i_k[(i - 1) * (n_sample - 1) + k] < 1e-7 && W_i_k[(i - 1) * (n_sample - 1) + k] > -1e-7)
            {
                if (k % 2 == 0)
                    break;
            }
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

void weigth_grid_i(double **weight_grid, long double *wik, int col, int n_sample, double Hj)
{
    long double *vkj = element_time(n_sample, Hj);
    for (int j = 0; j < n_sample - 1; j++)
    {
        for (int k = 0; k < n_sample - 1; k++)
        {
            weight_grid[j][col] += vkj[k] * wik[j * (n_sample - 1) + k];
            // if(wik[j * (n_sample - 1) + k] < 1e-8)
            // break;
        }
    }
    free(vkj);
}

double **cumulatve_weight_v2(int n_sample, int grid_size, double upper_bound, double lower_bound, char *file_name, int log)
{
    //double *H = scale_time(grid_size, upper_bound, lower_bound, file_name);
    double * H;
    if(log) 
        H = generate_logarithmic_scale(grid_size, upper_bound, lower_bound, file_name); 
    else
        H = generate_linear_scale(grid_size, upper_bound, lower_bound, file_name);
    long double *wik = Wik(n_sample);                                 // Calculate P_i,k values
    double **weight_grid = malloc((n_sample - 1) * sizeof(double *)); // Allocate memory for the weight grid
    for (int i = 0; i < n_sample - 1; i++)
        weight_grid[i] = calloc(grid_size + 2, sizeof(double)); // Initialize each row of the weight grid
    for (int i = 1; i < grid_size + 1; i++)
        weigth_grid_i(weight_grid, wik, i, n_sample, H[i - 1]);
    // Free the memory allocated for P_i,k values
    weigth_grid_i(weight_grid, wik, grid_size + 1, n_sample, INFINITY);
    free(H);
    free(wik);
    return weight_grid;
}