#ifndef SFS_H
#define SFS_H

typedef struct sfs
{
    double *training;
    double *test;
    int sfs_length;
    int n_haplotypes;
    double training_size;
    int singleton;
    int troncation;
    int oriented;
} SFS; 


SFS int_sfs(char *sfs_file, int oriented, int troncation, int singleton, int num_blocks);
void clear_sfs(SFS sfs);

#endif // SFS_H