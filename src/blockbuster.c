#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <lapacke.h>
#include <stddef.h>
#include <cblas.h>
#include <omp.h>
#include "blockbuster.h"
#include "blockbuster_grid.h"
#include "linear_regression.h"

/**
 * Generates the next combination of breakpoints in lexicographic order within a specified grid.
 *
 * @param brk            Array representing the current combination of breakpoints. breakpoints are index of times of changes from which cumulative branch lengthes have been computed
 *                       This array will be modified in place to hold the next combination.
 * @param nb_breakpoints The number of breakpoints (i.e., the size of the `brk` array).
 * @param grid_size      The size of the grid, representing the maximum index range for breakpoints.
 *
 * @return               1 if there are more combinations to generate,
 *                       0 if the last combination has been reached.
 */
int recursive_bk_combination(int *brk, int nb_breakpoints, int grid_size)
{
    // Start with the last breakpoint in the array
    int i = nb_breakpoints - 1;
    // Attempt to increment the last breakpoint by 1
    brk[i] += GRIDREFINE;
    // Adjust previous breakpoints if the current one exceeds allowable range
    while (i >= 0 && brk[i] / GRIDREFINE > grid_size + 1 - nb_breakpoints + i)
    {
        // Move to the previous breakpoint in the array
        i--;
        // Increment the previous breakpoint
        brk[i] += GRIDREFINE;
    }
    // For any indices after `i`, reset them to maintain increasing order
    for (int j = i + 1; j < nb_breakpoints; j++)
        brk[j] = brk[j - 1] + GRIDREFINE;
    // If the first breakpoint reaches its upper limit, all combinations have been generated
    if (i == 0 && brk[0] / GRIDREFINE == grid_size + 1 - nb_breakpoints)
        return 0;
    // Return 1 to indicate more combinations are available
    return 1;
}


/**
 * Generates and evaluates all possible breakpoint combinations to find the optimal solution
 * with the highest log-likelihood, saving each intermediate solution and returning the best one.
 *
 * @param nb_breakpoints  Number of breakpoints (population size changes) to consider.
 * @param sfs_length      Length of the Site Frequency Spectrum (SFS).
 * @param cumul_weight    Cumulative weights matrix for branch lengths across time intervals.
 * @param sfs             Observed SFS data (theoretical and observed SFS).
 * @param grid_size       Total number of grid points available for placing breakpoints.
 * @param n_sample        Sample identifier or index used for saving solutions.
 *
 * @return                The `solution` structure containing the optimal combination of breakpoints
 *                        with the highest log-likelihood value.
 */
Solution generate_brk_combinations(int nb_breakpoints, SFS sfs, Time_gride tg)
{
    // Initialize the solution with specified number of breakpoints
    Solution sol = init_solution_size(nb_breakpoints);
    // Set the last breakpoint to the grid limit (end boundary)
    sol.breakpoints[nb_breakpoints] = GRIDREFINE * tg.grid_size + 1;
    // Flag for stopping the combination generation
    int arret = 1;
    // Calculate the initial solution with the given breakpoints
    system_resolution(&sol, sfs, tg);
    // Initialize a temporary solution to keep track of the best found solution
    Solution tmp_sol = copy_solution(sol);
    // Loop through all breakpoint combinations until `arret` is set to 0
    while (arret && nb_breakpoints > 0)
    {
        // Generate the next combination of breakpoints
        arret = recursive_bk_combination(sol.breakpoints, sol.nb_breakpoints, tg.grid_size);
        // Calculate the solution (log-likelihood and distance) for the new combination
        system_resolution(&sol, sfs, tg);
        // If the new solution has a higher log-likelihood and all theta values are positive, keep its
        if (isnan(tmp_sol.log_likelihood) || sol.log_likelihood > tmp_sol.log_likelihood)
        {
            clear_solution(tmp_sol);
            tmp_sol = copy_solution(sol);
        }
    }
    clear_solution(sol);
    // thetas_se(&tmp_sol, sfs_length, cumul_weight);
    // residues(&tmp_sol, cumul_weight, sfs_length, sfs[0]);
    // Return the best solution with the optimal breakpoint combination for a givent nulber of changes in population size
    return tmp_sol;
}

// int check_new(solution sol, solution sol_initiale, int breakpoint)
// {
//     if(sol.breakpoints[breakpoint] > sol.breakpoints[breakpoint + 1] || sol.breakpoints[breakpoint] < sol.breakpoints[breakpoint - 1])
//         return 0;
//     return 1;
// }

Solution refine_solution_b(Solution sol_initiale, int b, SFS sfs, Time_gride tg, int sign)
{
    Solution sol1 = copy_solution(sol_initiale);
    Solution solm = copy_solution(sol_initiale);
    sol1.breakpoints[b]  += sign;
    system_resolution(&sol1, sfs, tg);
    
    while(sol1.log_likelihood > solm.log_likelihood){
        clear_solution(solm);
        solm = copy_solution(sol1);
        sol1.breakpoints[b]  += sign;
        system_resolution(&sol1, sfs, tg);
    }
    clear_solution(sol1);
    // clear_solution(sol2);
    return solm; 
}


void refine_solution(Solution *sol_initiale, SFS sfs, Time_gride tg)
{
    for(int i =0; i < 100; i ++)
    {   
        for(int b = sol_initiale->nb_breakpoints - 1; b >= 0; b --)
        {
            Solution solp = refine_solution_b(*sol_initiale, b, sfs, tg, -1);
            Solution solm = refine_solution_b(*sol_initiale, b, sfs, tg, +1);
            if(sol_initiale->log_likelihood < solm.log_likelihood)
            {
                clear_solution(*sol_initiale);
                *sol_initiale = copy_solution(solm);
            }
            else{
                if(sol_initiale-> log_likelihood < solp.log_likelihood)
                {
                    clear_solution(*sol_initiale);
                    *sol_initiale = copy_solution(solp);
                }
            }
            clear_solution(solp);
            clear_solution(solm);
        }
    }
}

/**
 * Finds the best scenarios by generating combinations of breakpoints and solving
 * the regression for each configuration to estimate population mutation rates (thetas).
 *
 * @param sfs_length   Length of the observed Site Frequency Spectrum (SFS).
 * @param cumul_weight Pointer to a 2D array of cumulative weights for regression.
 * @param sfs          Pointer to a 2D array where:
 *                      - sfs[0] is the training SFS used for regression.
 *                      - sfs[1] is the test SFS used for validation.
 * @param grid_size    size of the discrete gride of time points
 * @param n_sample     sample size
 * @param changes      Maximum number of allowed breakpoint changes.
 *
 * @return             Array of solutions containing the best solution for each number of 
 *                     population size changes, from 0 to `changes`. Each solution includes
 *                     estimated thetas and associated metrics (log-likelihood, distance).
 *
 * Steps:
 * 1. Calculate `const_ren`, the proportion of sites represented in the training SFS (`sfs[0]`).
 *    This is used to adjust scaling if needed.
 * 3. Iteratively generate scenarios for `nb_breakpoints` from 0 to `changes`:
 *    - Use `generate_brk_combinations` to compute the best solution for each number of breakpoints.
 * 4. Return the array of solutions for all configurations.
 */
Solution *find_scenario(SFS sfs, Time_gride tg, int changes)
{
    // Allocate memory for storing solutions for all configurations of breakpoints
    Solution *liste_solution = malloc(sizeof(Solution) * (changes + 1));

    // Iterate through the number of breakpoints from 0 to `changes`
    int nb_breakpoints = 0;
    while (nb_breakpoints <= changes)
    {
        // Generate the best solution for the current number of breakpoints
        liste_solution[nb_breakpoints] = generate_brk_combinations(nb_breakpoints, sfs, tg);
        if (nb_breakpoints >= 1)
        {
            refine_solution(&liste_solution[nb_breakpoints], sfs, tg);
        }
        printf("%d ", nb_breakpoints);
        nb_breakpoints++;
    }
    // thetas_se(&liste_solution[0], sfs_length, cumul_weight);
    // residues(&liste_solution[0], cumul_weight, sfs_length, sfs[0]);
    // Return the list of solutions
    return liste_solution;
}


// void copy_and_insert_in_sorted_array(int *breakpoints, int nb_breakpoint, int new_breakpoint, int *new_breakpoints)
// {
//     int i, j;
//     int inserted = 0;

//     // Parcours du tableau original et insertion dans le nouveau tableau
//     for (i = 0, j = 0; i < nb_breakpoint; i++, j++)
//     {
//         // Insérer le nouvel élément à la bonne position
//         if (!inserted && breakpoints[i] > new_breakpoint)
//         {
//             new_breakpoints[j] = new_breakpoint;
//             inserted = 1;
//             j++; // Pour éviter d'écraser l'élément suivant
//         }
//         // Copier l'élément du tableau original
//         new_breakpoints[j] = breakpoints[i];
//     }
// }


