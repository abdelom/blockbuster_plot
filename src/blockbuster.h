#ifndef PAST_FINDER_H
#define PAST_FINDER_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "linear_regression.h"

// Fonction pour trouver le scénario
Solution *find_scenario(SFS sfs, Time_gride tg, int changes);
// solution* recent_infrence(solution *list_solution, int changes, double **sfs, double ** cumul_weight, int sfs_length, int n_sample);
// void save_solution(solution sol, int n_sample, char *out_file, double const_ren, int sfs_length);

// // Fonction pour libérer la mémoire allouée pour l'intégrale des poids
// // void free_integral_grid(double **cumul_weight, int n_sample);
// solution init_solution_size(int nb_breakpoints);
// void system_resolution(solution *sol, double **sfs, double **cumul_weight, int sfs_length);
// // Fonction pour libérer la mémoire allouée pour la solution
// void clear_solution(solution sol);

#endif // PAST_FINDER_H
