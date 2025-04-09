#ifndef MSD_H
#define MSD_H

#include <stdint.h>

// Function declarations

double compute_MSD(double* coord1, double* coord2, int num_atoms);
double* compute_time_averaged_msd(const double* coord, const int64_t* timesteps, const int64_t num_timesteps, const int64_t deltaTimestep);
#endif  // MSD_HI
