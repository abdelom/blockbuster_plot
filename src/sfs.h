#ifndef SFS_H
#define SFS_H

/**
 * @brief Structure to store a Site Frequency Spectrum (SFS) and its related information.
 */
typedef struct sfs
{
    double *training;       // Array containing the training SFS.
    double *test;           // Array containing the test SFS (split from training).
    int sfs_length;         // Number of bins in the SFS.
    int n_haplotypes;       // Number of haplotypes (sample size = sfs_length + 1).
    double training_size;   // Proportion of data used for training (1 - 1/num_blocks).
    int singleton;          // Whether singletons are included (1 = yes, 0 = ignored).
    int troncation;         // Maximum number of bins to keep if truncation is applied.
    int oriented;           // Whether the SFS is oriented (1 = yes, 0 = folded/unoriented).
    double delta_time;    // Whether to use delta time splitting (1 = yes, 0 = no).
} SFS;

double generate_normal_random(double mu, double sigma);
SFS int_sfs(char *sfs_file, int oriented, int troncation, int singleton, int num_blocks, int delta_time);
SFS noise_sfs(SFS sfs);
void clear_sfs(SFS sfs);

#endif // SFS_H