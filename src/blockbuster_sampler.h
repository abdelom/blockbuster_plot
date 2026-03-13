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

#ifndef bls_H
#define bls_H

#include "sfs.h"
#include "blockbuster_grid.h"
#include "linear_regression.h"
#include "blockbuster.h"

Solution *bootstrap(SFS sfs, Time_gride tg, int epochs, int repeats);
void format_bootstrap(Solution *list_solution, int repeats, double mutation_rate, double genome_length, double generation_time, const char *prefix);

#endif