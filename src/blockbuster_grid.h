#ifndef PAST_GRID_H // Si ce fichier n'a pas encore été inclus
#define PAST_GRID_H // alors définir PAST_GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "sfs.h"
#include <omp.h>


typedef struct Time_gride
{
    double *time_scale;
    double **cumulative_bl;
    int grid_size;
} Time_gride; 

#define GRIDREFINE 1000

Time_gride init_time_grid(SFS sfs, int grid_size, double ub, double  lb, char *outfile);
Time_gride init_time_grid_H(int n_haplotypes, int grid_size, double *H);
void save_cumulated_weight(int sfs_length, int grid_size, double **matrix, char *filename);
void clear_time_grid(Time_gride tg, int sfs_length);
double* generate_linear_scale(int grid_size, double upper_bound, double lower_bound, char *file_name);
double* generate_logarithmic_scale(int grid_size, double upper_bound, double lower_bound, char *file_name);
#endif // Fin de la condition PAST_GRID