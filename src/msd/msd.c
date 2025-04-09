#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

double* compute_time_averaged_msd(const double* coordinates, const int64_t* timesteps, const int64_t num_timesteps, const int64_t timestep_difference) {
    // Number of windows must be computed with the formula in case timesteps are missing.
    // Having holes does not decrease the number of windows
    // Starting late does decrease the number of windows!
    // The maximum length possible is the last timestep divided by timestep_difference rounded up
    int64_t total_n_windows = (timesteps[num_timesteps-1] - timesteps[0] + timestep_difference - 1 ) / timestep_difference;
    total_n_windows+=1; // Zero length window. I expect it to be zero and compute that as a check

    double* ave_msd = (double*) calloc(total_n_windows, sizeof(double));
    if (!ave_msd) {
        fprintf(stderr, "Failed to allocate memory for ave_msd\n");
        exit(-1);
    }
    int64_t* window_counters = (int64_t*) calloc(total_n_windows, sizeof(int64_t));
    if (!window_counters) {
        fprintf(stderr, "Failed to allocate memory for the number of windows\n");
        exit(-1);
    }

    // `t` is the time index
    // `w_deltaT` is the length of the window
    for(int64_t t_start=0; t_start<num_timesteps; t_start++) {
        #pragma omp parallel for
        for(int64_t t_end=t_start+1; t_end<num_timesteps; t_end++) {
            int64_t time_difference = timesteps[t_end]-timesteps[t_start];
            int64_t idx_deltaTime = time_difference / timestep_difference; // [Assumption] Timesteps are all multiples of timestep_difference

            // Let's compute the MSD for the single particle
            const double dx = coordinates[3*t_end]   - coordinates[3*t_start];
            const double dy = coordinates[3*t_end+1] - coordinates[3*t_start+1];
            const double dz = coordinates[3*t_end+2] - coordinates[3*t_start+2];
            const double msd = dx*dx + dy*dy + dz*dz;
            ave_msd[idx_deltaTime] += msd;
            window_counters[idx_deltaTime] +=1;
        }
        //break;
    }

    // Normalize by the number of windows computed (time average)
    for(int64_t t=0; t<total_n_windows; t++) {
        if (window_counters[t]!=0) {
            ave_msd[t] /= (double) window_counters[t];
        } //else if(t>0) { printf("No match for window size=%ld\n", t); }
    }
    free(window_counters);
    return ave_msd;
}


double compute_MSD(const double* coordinates1, const double* coordinates2, const int num_atoms) {
    // Computes the MSD averaged over all particles at a given timestep.
    // It assumes a linear array of coordinates with x, y, z for each atom.
    //  coordinates1: coordinates at time t for all atoms
    //  coordinates2: coordinates at time t+dT for all atoms
    //  num_atoms: number of atoms (the routine will assume 3*num_atoms coordinates)
    //  dT: time difference between coordinates1 and coordinates2
    double msd = 0.0;

    for(int i = 0; i < num_atoms; i++) {
        const double dx = coordinates1[3*i]   - coordinates2[3*i];
        const double dy = coordinates1[3*i+1] - coordinates2[3*i+1];
        const double dz = coordinates1[3*i+2] - coordinates2[3*i+2];
        msd += dx*dx + dy*dy + dz*dz;
    }
    return msd / ( (double) num_atoms );
}