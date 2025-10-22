#include "sfs.h"
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
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
 * Example usage:
 *     int size;
 *     double *sfs = readSFSFromFile("sfs_data.txt", &size);
 *     // Use `sfs` array here
 *     free(sfs);  // Free the allocated memory when done
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

 void fold_sfs(SFS *sfs)
{
    int n = sfs->sfs_length;

    // Combiner la première moitié avec la deuxième
    for (int i = 0; i < n / 2; i++) {
        sfs->training[i] += sfs->training[n - 1 - i];
    }

    // Nouvelle taille (moitié + reste si n est impair)
    int new_size = n / 2 + n % 2;

    // Réallouer training (et test si nécessaire)
    sfs->training = realloc(sfs->training, new_size * sizeof(double));
    // sfs->test     = realloc(sfs->test,     new_size * sizeof(double));

    sfs->sfs_length = new_size;
}


void replace_with_0(SFS *sfs)
{
    if(!sfs->singleton)
        sfs->training[0] = 0.;
    if(sfs->troncation && sfs->troncation < sfs->sfs_length && sfs->troncation > 5)
    {
        for(int i = sfs->troncation; i < sfs->sfs_length; i ++)
            sfs->training[i] = 0.;
    }
}

// void fold_sfs(SFS sfs, int grid_size)
// {
//     for (int i = 0; i < sfs.sfs_length / 2; i++)
//     {
//         // Combine the corresponding elements from the first and second halves of the SFS
//         sfs.training[i] += sfs.training[sfs.sfs_length - 1 - i];

//         // Adjust cumulative weights accordingly
//         // for (int j = 0; j < grid_size + 2; j++)
//         //     cumulative_weight[i][j] += cumulative_weight[sfs_length - 1 - i][j];
//     }
//     sfs.sfs_length = sfs.sfs_length / 2 + sfs.sfs_length % 2;
//     // cumulative_weight = realloc(cumulative_weight, sfs_length/2 + sfs_length%2);
// }


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
 void singleton_erased(SFS *sfs)
{
    int n = sfs->sfs_length;

    // Décaler les valeurs (effacer le singleton)
    if(sfs->training[0] == 0)
        {
            for (int i = 0; i < n - 1; i++)
            sfs->training[i] = sfs->training[i + 1];

            // Réallocation avec la nouvelle taille
            sfs->training = realloc(sfs->training, (n - 1) * sizeof(double));
            // Mettre à jour la longueur
            sfs->sfs_length = n - 1;
        }
}

// void singleton_erased(SFS sfs)
// {
//     int n = sfs.sfs_length;

//     for (int i = 0; i < n - 1; i++)
//         *(sfs.training)[i] = *(sfs.training)[i + 1];
//     //     for (int j = 0; j < grid_size + 2; j++)
//     //         (*cumulative_weight)[i][j] = (*cumulative_weight)[i + 1][j];
//     // }

//     // Ajustement des tailles
//     // *cumulative_weight = realloc(*cumulative_weight, (n - 1) * sizeof(double *));
//     *(sfs.training) = realloc((*(sfs.training), (n - 1) * sizeof(double));
//     *(sfs.test) = realloc((*(sfs.test), (n - 1) * sizeof(double));

//     sfs.sfs_length = n - 1;
// }

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

SFS int_sfs(char *sfs_file, int oriented, int troncation, int singleton, int num_blocks)
{
    SFS sfs;
    sfs.sfs_length = 0;
    sfs.troncation = troncation;
    sfs.oriented = oriented;
    sfs.singleton = singleton;
    // Lecture du fichier et application des opérations sur SFS
    sfs.training = readSFSFromFile(sfs_file, &(sfs.sfs_length));
    sfs.n_haplotypes = sfs.sfs_length + 1; // n_sample doit être spécifié
    replace_with_0(&sfs);
    if(!oriented)
        fold_sfs(&sfs);
    if(troncation && troncation <  sfs.sfs_length && troncation > 5){
        // sfs_troncation(sfs, cumul_weight, size, args.grid_size * GRIDREFINE, args.troncation);
        sfs.sfs_length = troncation;
        sfs.training = realloc(sfs.training, sfs.sfs_length * sizeof(double));
    }
    if(!singleton)
        singleton_erased(&sfs);
    sfs.training_size = 1. - 1. / (double) num_blocks;
    if(num_blocks == 1)
        sfs.training_size = 1.;
    sfs.test = test_split(sfs, num_blocks);
    // for (int i = 0; i < sfs.sfs_length; i++)
    //     printf("%f %f \n", sfs.training[i], sfs.test[i]);
    return sfs;
}

void clear_sfs(SFS sfs)
{
    free(sfs.test);
    free(sfs.training);
}