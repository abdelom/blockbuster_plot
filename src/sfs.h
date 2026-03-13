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


#ifndef SFS_H
#define SFS_H

/**
 * @brief Structure to store a Site Frequency Spectrum (SFS) and its related information.
 */
typedef struct sfs
{
    double *training;       // Array containing the training SFS.
    double *test;           // Array containing the test SFS (split from training).
    int sfs_length;         // Number of bins in the SFS.
    int n_haplotypes;       // Number of haplotypes (sample size = sfs_length + 1).
    double training_size;   // Proportion of data used for training (1 - 1/num_blocks).
    int singleton;          // Whether singletons are included (1 = yes, 0 = ignored).
    int troncation;         // Maximum number of bins to keep if truncation is applied.
    int oriented;           // Whether the SFS is oriented (1 = yes, 0 = folded/unoriented).
    double delta_time;    // Whether to use delta time splitting (1 = yes, 0 = no).
} SFS;

double generate_normal_random(double mu, double sigma);
SFS int_sfs(char *sfs_file, int oriented, int troncation, int singleton, int num_blocks, int delta_time);
SFS noise_sfs(SFS sfs);
void clear_sfs(SFS sfs);

#endif // SFS_H