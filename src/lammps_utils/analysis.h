//
// Created by gu on 27/03/25.
//
#ifndef LAMMPSDATAANALYSIS_H
#define LAMMPSDATAANALYSIS_H

#include <stdint.h>

void compute_CoM(const double* frame, const int64_t num_atoms, double* com_coord);

#endif //LAMMPSDATAANALYSIS_H
