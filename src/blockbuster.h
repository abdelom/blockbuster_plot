/*
        Copyright (C) 2026 A Omarjee

        This program is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public License
        as published by the Free Software Foundation; either version 2.1
        of the License, or (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
        for more information, please contact Abdelmajid Omarjee <abdelmajid.omarjee@mnhn.fr>
*/

#ifndef PAST_FINDER_H
#define PAST_FINDER_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "linear_regression.h"
#include "linear_regression_f.h"

// Fonction pour trouver le scénario
Solution *find_scenario(SFS sfs, Time_gride tg, int changes, int r);
// Solution generate_brk_combinations_f(int nb_breakpoints, SFS sfs, Time_gride tg, Flag f);
Solution *find_scenario_f(SFS sfs, Time_gride tg, Flag flag);
Solution generate_brk_combinations(int nb_breakpoints, SFS sfs, Time_gride tg, int r);
void refine_solution(Solution *sol_initiale, SFS sfs, Time_gride tg, int r);
int logratio_cumulative_test(const Solution *solutions, int k, double alpha);
void print_solution(Solution sol, Time_gride tg);

// solution* recent_infrence(solution *list_solution, int changes, double **sfs, double ** cumul_weight, int sfs_length, int n_sample);
// void save_solution(solution sol, int n_sample, char *out_file, double const_ren, int sfs_length);

// // Fonction pour libérer la mémoire allouée pour l'intégrale des poids
// // void free_integral_grid(double **cumul_weight, int n_sample);
// solution init_solution_size(int nb_breakpoints);
// void system_resolution(solution *sol, double **sfs, double **cumul_weight, int sfs_length);
// // Fonction pour libérer la mémoire allouée pour la solution
// void clear_solution(solution sol);

#endif // PAST_FINDER_H
