//
// Created by Guglielmo Grillo on 28/01/25.
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


#include "IO.h"
#include "integrator.h"
#include "physics.h"


void loadParamsFromFile(char *fileName, drivingFParams* dParams, IIOrdDiffEq* odeParams, integrator* integrator) {
    char key[10]; // The example shows key made of maximum 2 characters. I choose 10 just in case.
    double value;

    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("Error opening file %s\n", fileName);
        exit(1);
    }

    double k =  NAN,
           G =  NAN,
           F =  NAN,
           W =  NAN,
           x0 = NAN,
           v0 = NAN,
           dt = NAN,
           T  = NAN; // Declare all the variable to read.
    const double m = 1; // This is known from the text problem. I'll hardcode it for completeness.

    // Keep reading until find less than two fields (e.g. if file is finished)
    while ( fscanf(fp, "%s %lf", key, &value) == 2) {
        if (strcmp(key, "k") == 0) {
            k = value;
        } else if (strcmp(key, "G") == 0) {
            G = value;
        } else if (strcmp(key, "F") == 0) {
            F = value;
        } else if (strcmp(key, "W") == 0) {
            W = value;
        } else if (strcmp(key, "x0") == 0) {
            x0 = value;
        } else if (strcmp(key, "v0") == 0) {
            v0 = value;
        } else if (strcmp(key, "dt") == 0) {
            dt = value;
        } else if (strcmp(key, "T") == 0) {
            T = value;
        } else {
            printf("Unknown key-value: %s %lf\n", key, value);
            exit(1);
        }
    }

    bool missing_values = false;
    if(isnan(k)) {
        missing_values = true;
        printf("k is missing\n");
    }
    if(isnan(G)) {
        missing_values = true;
        printf("G is missing\n");
    }
    if(isnan(F)) {
        missing_values = true;
        printf("F is missing\n");
    }
    if( isnan(W)) {
        missing_values = true;
        printf("W is missing\n");
    }
    if( isnan(x0)) {
        missing_values = true;
        printf("x0 is missing\n");
    }
    if( isnan(v0)) {
        missing_values = true;
        printf("v0 is missing\n");
    }
    if( isnan(dt)) {
        missing_values = true;
        printf("dt is missing\n");
    }
    if( isnan(T)) {
        missing_values = true;
        printf("T is missing\n");
    }
    if(missing_values) {
        printf("Missing values! The program will quit.\n");
        exit(EXIT_FAILURE);
    }

    dParams->F = F;
    dParams->W = W;

    odeParams->a = -k/m;
    odeParams->b = -G/m;
    odeParams->t = 0.0;
    odeParams->x = x0;
    odeParams->v = v0;

    integrator->T = T;
    integrator->dt = dt;

    return;
}

void print_status(IIOrdDiffEq* ode) {
    printf("%lf, %lf, %lf, %lf\n", ode->t, ode->x, ode->v, idealHOenergy(ode));
}