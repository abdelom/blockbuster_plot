#ifndef PAST_GRID_H // Si ce fichier n'a pas encore été inclus
#define PAST_GRID_H // alors définir PAST_GRID_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
// #include <mpfr.h>
#include <omp.h>

//mpfr_t *coalescents_rate_hp(int n_sample);
//mpfr_t **init_mpfrt_matrix(int n_sample);
// void ratio2(mpfr_t rat, mpfr_t a, mpfr_t b, mpfr_t c, int precision);
// mpfr_t **hypoexponential_consts_hp(mpfr_t *rate, int n_sample);
// void exp_rate_time(mpfr_t result, mpfr_t rate, double h);
// double integral_calcul(mpfr_t *constants, mpfr_t *rate, int k, double hj, int precision);
// void free_mpfr_arrays(mpfr_t **constants, mpfr_t *rate, int n_sample);
// double **integral_grid_calcul(int n_sample, int *grid_size, double upper_bound, double lower_bound, char *file_name);
// double *k_pik_2(int n);
void matrix_multiply(int n, int t, double *A, double **B, double **C);
void save_cumulated_weight(int n_sample, int grid_size, double **matrix, const char *filename);
void free_integral_grid(double **cjk_grid, int n_sample);
double* generate_logarithmic_scale(int grid_size, double upper_bound, double lower_bound, char *file_name);
double* generate_linear_scale(int grid_size, double upper_bound, double lower_bound, char *file_name);
//double **cumulatve_weight(int n_sample, int *grid_size, double upper_bound,  double lower_bound, char *file_name);
double **cumulatve_weight_v2(int n_sample, int grid_size, double *H);
//double *coalescents_rate(int n_sample);
//double **hypoexponential_consts(double *rate, int n_sample);
long double *Wik(int n_sample);
#endif // Fin de la condition PAST_GRID