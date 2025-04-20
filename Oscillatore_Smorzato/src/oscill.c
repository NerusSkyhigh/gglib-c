#include <stdlib.h>
#include <stdio.h>

#include "integrator.h"
#include "IO.h"
#include "physics.h"


int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }
    char *filename = argv[1];


    drivingFParams* dFP = malloc(sizeof(drivingFParams));
    IIOrdDiffEq* ode    = malloc(sizeof(IIOrdDiffEq));
    integrator* RKI     = malloc(sizeof(integrator)); // Runge-Kutta Integrator

    loadParamsFromFile(filename, dFP, ode, RKI);

    ode->dParams = dFP;
    ode->f = &sinusoidal;
    RKI->integrate = &runge_kutta_integrator;

    /*
    printf("DRIVING_PARAMS\n");
    printf("\tdFP->F=%lf \n\tdFP->W=%lf\n", dFP->F, dFP->W);

    printf("ODE PARAMS:\n");
    printf("\tode->a=%lf \n\tode->b=%lf \n\tode->t=%lf \n\tode->x=%lf \n\tode->v=%lf\n", ode->a, ode->b, ode->t, ode->x, ode->v);

    printf("RKI PARAMS\n");
    printf("\tRKI->dt=%lf \n\tRKI->T=%lf\n", RKI->dt, RKI->T);
    */

    while( ode->t <= RKI->T ) {
        RKI->integrate(ode, RKI);
        ode->t += RKI->dt;
        print_status(ode);

    }


    free(dFP);
    free(ode);
    free(RKI);

    return 0;
}
