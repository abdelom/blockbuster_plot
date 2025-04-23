typedef struct
{
    int n_samples;
    float l_genome;
    double theta;
    int oriented;
    int seed;
    char outfile[100];
    double *parameters;
    int n_parameters;
    int noised;
    int header;
} Args;


double *sfs_infinite(Args args, double **cumulated_branch_lengthes);
double *folded_sfs(double *sim_sfs, int sfs_length);