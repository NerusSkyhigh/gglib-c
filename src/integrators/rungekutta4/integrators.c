/**
 * @author Guglielmo Grillo
 */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#include "integrators.h"

/**
 * @brief Runge Kutta integrator of order 4
 * @param ps PhysicsSystem which contains the state of the system in a specific frame
 * @param integrator parameters of the integrator
 */
void runge_kutta_integrator(PhysicsSystem* ps, const Integrator* integrator) {
    const double dt = integrator->dt;
    const double hdt = 0.5 * dt; // Half dt; (need advice on this)

    // Store previous positions and velocities in the data buffer
    RK4_data* data = (RK4_data*) integrator->_data;
    memcpy(data->bx, ps->x, 3*ps->N*sizeof(double));
    memcpy(data->bv, ps->v, 3*ps->N*sizeof(double));

    // Pointers to the buffer to simplify notation
    double* xk1=data->xk1;    double* vk1=data->vk1;
    double* xk2=data->xk2;    double* vk2=data->vk2;
    double* xk3=data->xk3;    double* vk3=data->vk3;
    double* xk4=data->xk4;    double* vk4=data->vk4;

    // NOTE: As the function inside integrator takes only a PhysicsSystem
    // as a parameter, t, x, and v need to be update inside the struct.
    // The original value will be stored in the data_ field
    data->bt = ps->t;
    double* x = data->bx;
    double* v = data->bv;


    // We are solving two different coupled ode at the same time.
    // I assume that solving them independently is ok.
    //for (int64_t p=0; p<3*ps->N; p++) {
    //    xk1[p] = v[p];
    //    vk1[p] = integrator->f(ps);
    //}
    memcpy(xk1, v, 3*ps->N*sizeof(double));
    integrator->f(ps, vk1);

    ps->t = data->bt+hdt;
    for (int64_t p=0; p<3*ps->N; p++) {
        //xk2[p] = v[p] + hdt*vk1[p];

        // Update the positions in the physics to the middle step
        ps->x[p] = x[p] + hdt*xk1[p];
        ps->v[p] = v[p] + hdt*vk1[p];
    }
    memcpy(xk2, ps->v, 3*ps->N*sizeof(double));
    integrator->f(ps, vk2);


    //ps->t = data->bt+hdt; // Already set from previous state
    for (int64_t p=0; p<3*ps->N; p++) {
        //xk3[p] = v[p] + hdt*vk2[p];
        ps->x[p] = x[p] + hdt*xk2[p];
        ps->v[p] = v[p] + hdt*vk2[p];
    }
    memcpy(xk3, ps->v, 3*ps->N*sizeof(double));
    integrator->f(ps, vk3);

    ps->t = data->bt+dt;
    for (int64_t p=0; p<3*ps->N; p++) {
        //xk4[p] = v[p] + dt*vk3[p];
        ps->x[p] = x[p] + dt*xk3[p];
        ps->v[p] = v[p] + dt*vk3[p];
    }
    memcpy(xk4, ps->v, 3*ps->N*sizeof(double));
    integrator->f(ps, vk4);

    // Update the original timestep in the PhysicsSystem struct
    // The consts help the compiler optimize the division by
    // compiling with the division already done and performing
    // a multiplication instead
    const double inv3 = 1.0/3.0;
    const double inv6 = 1.0/6.0;
    for (int64_t p=0; p<3*ps->N; p++) {
        //ps->x[p] = x[p] + dt* (xk1[p]/6. + xk2[p]/3. + xk3[p]/3. + xk4[p]/6.);
        //ps->v[p] = v[p] + dt* (vk1[p]/6. + vk2[p]/3. + vk3[p]/3. + vk4[p]/6.);
        ps->x[p] = x[p] + dt* ( (xk1[p]+xk4[p])*inv6 + (xk2[p]+xk3[p])*inv3 );
        ps->v[p] = v[p] + dt* ( (vk1[p]+vk4[p])*inv6 + (vk2[p]+vk3[p])*inv3 );
    }
    ps->t = data->bt+dt;


}


/**
 * @brief Init the integrator with a runge kutta integrator of 4th order
 * @param rk4i pointer to an Integrator with fields dt, t and f already fixed
 * @param ps pointer to a PhysicsSystem used to compute the content of _data
 */
void init_RK4(Integrator* rk4i, const PhysicsSystem* ps) {
    printf("Trying to add runge kutta\n");
    rk4i->integrate = runge_kutta_integrator;
    rk4i->_data = malloc(sizeof(RK4_data));
    RK4_data* data = (RK4_data*) rk4i->_data;

    printf("Malloc problem\n");
    // To enforce memory cohesion, a single linear array is allocated
    // It will contain xk and vk for each particle in three dimensions.
    // There are also 4 different xk and vk. Furthermore, we have to
    // store the partial values for position and velocity.
    // In total, we have 2*4*3*N+2*3N elements with the following layout:
    //  xk1_x1, xk1_y1, xk1_z1, xk1_x2, ...,
    //      vk1_x1, vk1_y1, vk1_z1, vk1_x2, ...,
    //          xk2_x1, xk2_y1, xk2_z1, xk2_x2, ...
    //              ...
    //                  ..., vk4_xN, vk4_yN, vk4_zN, ...
    //                          x1, y1, z1, ..., vxN, vyN, vzN
    // That is, there are 10 vectors of length 3N
    const int64_t length = 3*ps->N;
    data->xk1 = calloc(10*length, sizeof(double));
    data->vk1 = data->xk1 + length; // From 3N to 6N
    data->xk2 = data->vk1 + length; // From 6N to 9N
    data->vk2 = data->xk2 + length; // From 9N to 12N
    data->xk3 = data->vk2 + length; // From 12N to 15N
    data->vk3 = data->xk3 + length; // From 15N to 18N
    data->xk4 = data->vk3 + length; // From 18N to 21N
    data->vk4 = data->xk4 + length; // From 21N to 24N
    data->bx  = data->vk4 + length; // From 24N to 27N
    data->bv  = data->bx  + length; // From 27N to 30N
}

/**
 * @brief Free the rk4i struct and its content.
 * @param rk4i pointer to the struct to free
 */
void free_RK4(Integrator* rk4i) {
    const RK4_data* data = (RK4_data*) rk4i->_data;
    free(data->xk1); // As it is a linear array I have to free only the first pointer
    free(rk4i->_data);
}


