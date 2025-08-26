#ifndef LINEAR_REGRESSION_F_H
#define LINEAR_REGRESSION_F_H


typedef struct
{
    int n_theta;
    int n_theta_fixed;
    double * thetas;
    double * sfs;
} Flag;

Flag init_flag(int sfs_length, double *thetas , double n_theta);
double * init_system_f(SFS sfs, Flag *flag, System system);
void system_resolution_f(Solution *sol, SFS sfs, Time_gride tg, Flag flag);
#endif