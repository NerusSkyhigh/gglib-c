/**
 * @author Guglielmo Grillo
 * @brief This header contains the function necessary to perform the numerical integration of the equation of motion
 */
#pragma once

#include <stdint.h>

/**
 * @brief PhysicsSystem_ struct. Contains the state of the system in a specific frame
 * @param x pointer to a linear array of size 3N containing the coordinates of the system: [x1, y1, z1, x2, y2, z2, x3..., xN, yN, zN]
 * @param v pointer to a linear array of size 3N containing the coordinates of the system: [x1, y1, z1, x2, y2, z2, x3..., xN, yN, zN]
 * @param m pointer to a linear array of size N containing the masses of the particles
 * @param N number of particles in the system
 */
typedef struct PhysicsSystem_ {
    double t;
    double* x;
    double* v;
    double* m;
    int64_t N;
} PhysicsSystem;

typedef void (*acceleration_f)(const PhysicsSystem*, double*);
/**
 * @brief Integrator_ struct. Contains the data used by the integrator
 * @param dt timestep for each integrator step
 * @param t current time of the integrator
 * @param f pointer to a function that computes the derivative of the velocity and stores it in the second argument.
 * @param integrate pointer to the chosen inte grator which takes the current state of the system
 * @param data_ contains data used by the specific integrator chosen.
 */
typedef struct Integrator_ {
    double dt;
    acceleration_f f;
    void (*integrate)(PhysicsSystem*, const struct Integrator_*);
    void* _data;
} Integrator;
typedef void (*integrate_f)(PhysicsSystem*, const Integrator*);


/**
 * Runge Kutta 4th Order functions
 */

/**
 * @brief _data for the Rk4 integrator
 * @param xk1, xk2, xk3, xk4 coefficients for the position integration
 * @param vk1, vk2, vk3, vk4 coefficients for the velocity integration
 * @param bx, bv, bt buffers for of x, v, and t
 */
typedef struct RK4_data_ {
    double* xk1;    double* vk1;
    double* xk2;    double* vk2;
    double* xk3;    double* vk3;
    double* xk4;    double* vk4;
    double* bx;     double* bv;
    double bt;
} RK4_data;
void init_RK4(Integrator*, const PhysicsSystem*);
void free_RK4(Integrator*);
void runge_kutta_integrator(PhysicsSystem*, const Integrator*);