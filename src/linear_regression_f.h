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