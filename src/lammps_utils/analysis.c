//
// Created by gu on 27/03/25.
//

#include <stdint.h>
#include "parser.h"

void compute_CoM(const double* frame, const int64_t num_atoms, double* com_coord) {
    // Set to zero for accumulation
    com_coord[0] = 0.;
    com_coord[1] = 0.;
    com_coord[2] = 0.;

    for (int i = 0; i < num_atoms; i++) { // TODO: Transform this into a `size_t` object
        // Let's divide by number of atoms to avoid overflow
        com_coord[0] += frame[3*i]   / ( (double) num_atoms);
        com_coord[1] += frame[3*i+1] / ( (double) num_atoms);
        com_coord[2] += frame[3*i+2] / ( (double) num_atoms);
    }
    //printf("%lf %lf %lf\n", com_coord[0], com_coord[1], com_coord[2]);
}