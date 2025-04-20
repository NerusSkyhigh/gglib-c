//
// Created by Guglielmo Grillo on 28/01/25.
//
#ifndef IO_H
#define IO_H

#include "integrator.h"

void loadParamsFromFile(char *fileName, drivingFParams* dParams, IIOrdDiffEq* odeParams, integrator* integrator);

void print_status(IIOrdDiffEq* ode);

#endif //IO_H
