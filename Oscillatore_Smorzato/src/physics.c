//
// Created by Guglielmo Grillo on 28/01/25.
//
#include <math.h>

#include "physics.h"
#include "integrator.h"


#include <stdio.h>

double sinusoidal(double t, double x, drivingFParams* dParams) {
    return dParams->F * sin( dParams->W * t);
}

double idealHOenergy(IIOrdDiffEq* ode) {
  // We exploid the assumption m=1 to use ode->a = -k/m = -k
  return 0.5*( -1.* ode->a * ode->x*ode->x + ode->v*ode->v);
}