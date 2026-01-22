#include "sfs.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
 * Reads the Site Frequency Spectrum (SFS) data from a file and returns it as a dynamically allocated array.
 * The function reads all floating-point numbers from the file, stores them as doubles, and returns the
 * array and its size through a pointer argument.
 *
 * @param filename The path to the file containing the SFS data, where each number represents an SFS entry.
 * @param size Pointer to an integer where the function will store the number of entries in the SFS.
 *
 * @return A dynamically allocated array of doubles containing the SFS data from the file. The caller is
 *         responsible for freeing this memory.
 *
 */
double *readSFSFromFile(const char *filename, int *size)
{
    // Open the file in read mode
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1); // Exit if file cannot be opened
    }
    // Count the number of floating-point numbers in the file
    float number;
    int count = 0;
    while (fscanf(file, "%f", &number) == 1)
        count++;
    // Reset file pointer to the beginning of the file to read the values into the array
    fseek(file, 0, SEEK_SET);
    // Allocate memory for the SFS array based on the counted number of entries
    double *sfs = (double *)malloc(count * sizeof(double));
    if (sfs == NULL)
    {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1); // Exit if memory allocation fails
    }
    // Read each float from the file into the array as a double
    int i = 0;
    while (fscanf(file, "%lf", &sfs[i]) == 1)
        i++;
    *size = count; // Set the size to the number of entries in the array
    // Close the file and return the allocated array
    fclose(file);
    return sfs;
}


/**
 * @brief Folds a Site Frequency Spectrum (SFS) to account for unpolarized data.
 *
 * This function combines the first half of the SFS with the second half
 * (mirrored), effectively folding the spectrum. If the length is odd,
 * the middle bin is preserved. The training array is resized accordingly.
 *
 * @param sfs Pointer to the SFS structure to modify.
 */
void fold_sfs(SFS *sfs)
{
    int n = sfs->sfs_length;

    // Combine first half with mirrored second half
    for (int i = 0; i < n / 2; i++) {
        sfs->training[i] += sfs->training[n - 1 - i];
    }

    // Compute new size (half plus middle if n is odd)
    int new_size = n / 2 + n % 2;

    // Reallocate training array to the new size
    sfs->training = realloc(sfs->training, new_size * sizeof(double));
    // Optionally, realloc test array if needed:
    // sfs->test = realloc(sfs->test, new_size * sizeof(double));

    sfs->sfs_length = new_size;
}

/**
 * @brief Replaces specific bins in the training SFS with 0 according to preprocessing rules.
 *
 * This function sets the first bin (singleton) to 0 if singletons are to be ignored.
 * It also zeroes out bins beyond the truncation threshold if truncation is applied.
 *
 * @param sfs Pointer to the SFS structure to modify.
 */
void replace_with_0(SFS *sfs)
{
    // Set singleton bin to 0 if singletons are ignored
    if (!sfs->singleton)
        sfs->training[0] = 0.0;

    // Set bins beyond the truncation point to 0
    if (sfs->troncation && sfs->troncation < sfs->sfs_length && sfs->troncation > 5)
    {
        for (int i = sfs->troncation; i < sfs->sfs_length; i++)
            sfs->training[i] = 0.0;
    }
}


/**
 * @brief Removes the singleton bin from the training SFS if its value is 0.
 *
 * This function shifts all values in the training SFS to the left to remove
 * the first bin (singleton), reallocates the array to the new size,
 * and updates the sfs_length accordingly.
 *
 * @param sfs Pointer to the SFS structure to modify.
 */
void singleton_erased(SFS *sfs)
{
    int n = sfs->sfs_length;

    // Check if the first bin (singleton) is 0
    if (sfs->training[0] == 0)
    {
        // Shift all values to the left to remove the singleton
        for (int i = 0; i < n - 1; i++)
            sfs->training[i] = sfs->training[i + 1];

        // Reallocate memory for the reduced array
        sfs->training = realloc(sfs->training, (n - 1) * sizeof(double));

        // Update the SFS length
        sfs->sfs_length = n - 1;
    }
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
    // srand(time(NULL));
    // Tirer au hasard un tiers des indices sans remise
    for (int i = 0; i < taille; i++)
    {
        int randIndex = rand() % (n - i);        // Tirer un indice aléatoire dans les indices restants
        tirages[i] = indices[randIndex];         // Stocker l'indice tiré
        indices[randIndex] = indices[n - i - 1]; // Remplacer l'indice tiré par le dernier non tiré
    }
    // Trier la liste des tirages
    qsort(tirages, taille, sizeof(int), comparer);
    free(indices);
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

double *test_split(SFS sfs, int frac)
{
    double n_sites = 0;
    for (int i = 0; i < sfs.sfs_length; i++)
        n_sites += sfs.training[i];
    int taille_tirages = n_sites / frac;
    int *tirages = malloc(sizeof(int) * taille_tirages);
    double *cumulated_sfs = malloc(sizeof(double) * sfs.sfs_length);
    double *sfs_test = malloc(sizeof(double) * sfs.sfs_length);
    if (frac == 1)
    {
        for (int i = 0; i < sfs.sfs_length; i++)
            sfs_test[i] = sfs.training[i];
        free(tirages);
        free(cumulated_sfs);
        return sfs_test;
    }
    site_sampling(n_sites, tirages, taille_tirages);
    cumulate_sfs(sfs.training, sfs.sfs_length, cumulated_sfs);
    difference_count(tirages, taille_tirages, cumulated_sfs, sfs.sfs_length, sfs_test);
    mettre_a_jour_count(sfs_test, sfs.sfs_length);
    for (int i = 0; i < sfs.sfs_length; i++)
        sfs.training[i] -= sfs_test[i];
    free(tirages);
    free(cumulated_sfs);
    return sfs_test;
}


/**
 * @brief Initializes an SFS structure from a file and applies preprocessing steps.
 *
 * This function reads the SFS from the input file, applies optional folding,
 * truncation, singleton removal, and splits it into training and test sets.
 *
 * @param sfs_file      Path to the input SFS file.
 * @param oriented      1 if the SFS is oriented, 0 if it should be folded.
 * @param troncation    Maximum number of bins to keep (if truncation is applied).
 * @param singleton     1 to keep singletons, 0 to ignore them.
 * @param num_blocks    Number of blocks for splitting the SFS into training/test sets.
 *
 * @return An initialized SFS structure containing training and test arrays and metadata.
 */
SFS int_sfs(char *sfs_file, int oriented, int troncation, int singleton, int num_blocks, int delta_time)
{
    SFS sfs;

    // Initialize basic parameters
    sfs.sfs_length = 0;
    sfs.troncation = troncation;       // If truncation = k, only the first k bins are kept
    sfs.oriented = oriented;           // 1 = oriented, 0 = folded/unoriented
    sfs.singleton = singleton;         // 1 = singletons included, 0 = ignored

    // Read SFS from file
    sfs.training = readSFSFromFile(sfs_file, &(sfs.sfs_length));
    sfs.n_haplotypes = sfs.sfs_length + 1; // Sample size

    // Preprocessing steps
    replace_with_0(&sfs);             // Replace ignored bins with 0
    if (!oriented)
        fold_sfs(&sfs);               // Fold the SFS if unoriented

    // Apply truncation if requested
    if (troncation && troncation < sfs.sfs_length && troncation > 5)
    {
        sfs.sfs_length = troncation;
        sfs.training = realloc(sfs.training, sfs.sfs_length * sizeof(double));
    }

    // Remove singleton bin if requested
    if (!singleton)
        singleton_erased(&sfs);

    // Determine training size proportion
    sfs.training_size = 1.0 - 1.0 / (double) num_blocks;
    if (num_blocks == 1)
        sfs.training_size = 1.0;

    // Split SFS into training and test sets
    sfs.test = test_split(sfs, num_blocks);
    if(delta_time)
        sfs.delta_time = 1e-1;
    return sfs;
}

SFS copy_sfs(SFS sfs)
{
    SFS new_sfs;
    new_sfs.sfs_length = sfs.sfs_length;
    new_sfs.n_haplotypes = sfs.n_haplotypes;
    new_sfs.training_size = sfs.training_size;
    new_sfs.singleton = sfs.singleton;
    new_sfs.troncation = sfs.troncation;
    new_sfs.oriented = sfs.oriented;
    new_sfs.training = malloc(sizeof(double) * sfs.sfs_length);
    new_sfs.test = malloc(sizeof(double) * sfs.sfs_length);
    for (int i = 0; i < sfs.sfs_length; i++)
    {
        new_sfs.training[i] = sfs.training[i];
        new_sfs.test[i] = sfs.test[i];
    }
    return new_sfs;
}


SFS noise_sfs(SFS sfs)
{
    SFS new_sfs;
    new_sfs.sfs_length = sfs.sfs_length;
    new_sfs.n_haplotypes = sfs.n_haplotypes;
    new_sfs.training_size = sfs.training_size;
    new_sfs.singleton = sfs.singleton;
    new_sfs.troncation = sfs.troncation;
    new_sfs.oriented = sfs.oriented;
    new_sfs.training = malloc(sizeof(double) * sfs.sfs_length);
    new_sfs.test = malloc(sizeof(double) * sfs.sfs_length);
    for (int i = 0; i < sfs.sfs_length; i++)
    {
        new_sfs.training[i] = generate_normal_random(sfs.training[i], sfs.training[i]);
        new_sfs.test[i] = new_sfs.training[i]; //generate_normal_random(sfs.training[i], sfs.training[i]); //generate_normal_random(sfs.test[i], sfs.test[i]);
    }
    return new_sfs;
}

double generate_normal_random(double mu, double sigma)
{
    double u1 = rand() / (double)RAND_MAX; // Génère un nombre aléatoire entre 0 et 1
    double u2 = rand() / (double)RAND_MAX;
    // Génère un autre nombre aléatoire entre 0 et 1
    double z = sqrt(-2 * log(u1)) * cos(2 * 3.1415926535898 * u2); // Transformation de Box-Muller
    return mu + z * sqrt(sigma);                                   // réalisation du variable aléatoire suivant une loie normale de paramètres mu et sigma
}


// SFS noise_sfs(SFS sfs)
// {
//     SFS new_sfs;
//     new_sfs.sfs_length = sfs.sfs_length;
//     new_sfs.n_haplotypes = sfs.n_haplotypes;
//     new_sfs.training_size = sfs.training_size;
//     new_sfs.singleton = sfs.singleton;
//     new_sfs.troncation = sfs.troncation;
//     new_sfs.oriented = sfs.oriented;
//     new_sfs.training = malloc(sizeof(double) * sfs.sfs_length);
//     new_sfs.test = malloc(sizeof(double) * sfs.sfs_length);
//     for (int i = 0; i < sfs.sfs_length; i++)
//     {
//         new_sfs.training[i] = generate_normal_random(sfs.training[i], sfs.training[i]);
//         new_sfs.test[i] = generate_normal_random(sfs.training[i], 20 *sfs.training[i]); //generate_normal_random(sfs.test[i], sfs.test[i]);
//     }
//     return new_sfs;
// }

void clear_sfs(SFS sfs)
{
    free(sfs.test);
    free(sfs.training);
}