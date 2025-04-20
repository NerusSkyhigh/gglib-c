//
// Created by Guglielmo Grillo on 28/01/25.
//

#include "integrator.h"

double inline _compile_ode(IIOrdDiffEq* ode, double t, double x, double v) {
    return ode->a * x + ode->b * v + ode->f(t, x, ode->dParams);
}

void runge_kutta_integrator(IIOrdDiffEq* ode, integrator* iParams) {
    // I'll set variables to static to avoid reserving space
    // at each function call. I'm not sure this is efficient.
    // An other way could be to use a arrays.
    static double xk1, xk2, xk3, xk4;
    static double vk1, vk2, vk3, vk4;
    double hdt = 0.5 * iParams->dt; // Half dt; (need advice on this)

    // Notice that we are solving two different coupled ode at the same time.
    // I'll solve them independently
    xk1 = ode->v;
	vk1 = _compile_ode(ode, ode->t, ode->x, ode->v);

    xk2 = ode->v + hdt*vk1;
    vk2 = _compile_ode(ode, ode->t+hdt, ode->x + hdt*xk1, ode->v + hdt*vk1);

    xk3 = ode->v + hdt*vk2;
    vk3 = _compile_ode(ode, ode->t+hdt, ode->x + hdt*xk2, ode->v + hdt*vk2);

    xk4 = ode->v + iParams->dt*vk3;
    vk4 = _compile_ode(ode, ode->t+iParams->dt, ode->x + iParams->dt*xk3, ode->v + iParams->dt*vk3);

    ode->x = ode->x + iParams->dt * ( xk1/6. + xk2/3. + xk3/3. + xk4/6. );
    ode->v = ode->v + iParams->dt * ( vk1/6. + vk2/3. + vk3/3. + vk4/6. );
}