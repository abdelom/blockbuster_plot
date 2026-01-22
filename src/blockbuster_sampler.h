
#ifndef bls_H
#define bls_H

#include "sfs.h"
#include "blockbuster_grid.h"
#include "linear_regression.h"
#include "blockbuster.h"

Solution *bootstrap(SFS sfs, Time_gride tg, int epochs, int repeats);
void format_bootstrap(Solution *list_solution, int repeats, double mutation_rate, double genome_length, double generation_time, const char *prefix);

#endif