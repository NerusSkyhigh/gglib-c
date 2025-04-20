//
// Created by Guglielmo Grillo on 28/01/25.
//

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

// function RKintegrator

/**
 * Contains the parameters the harmonic driving force
 *      f(t, x, drivingFParams) = F sin(w t)
 */
typedef struct drivingFParams {
    double F;
    double W;
} drivingFParams;

/**
 * Contains the parameters of a second order
 * differential equation of the form
 *      x'' = - k/m x -gamma/m x' + F sin(w t)
 *          = a x + b x' + f(t, x, drivingFParams)
 */
typedef struct IIOrdDiffEq {
    double a, b;
    drivingFParams* dParams;
    double (*f)(double, double, struct drivingFParams*);
    double t;
    double x;
    double v;
} IIOrdDiffEq;

/**
 * Contains the parameters of the integrator and the
 * corresponding integrator function
 */
typedef struct integrator {
    double dt;
    double T;
    void (*integrate)(IIOrdDiffEq*, struct integrator*);
} integrator;

double _compile_ode(IIOrdDiffEq* ode, double t, double x, double v);
void runge_kutta_integrator(IIOrdDiffEq* ode, integrator* iParams);


#endif //INTEGRATOR_H
