//
// Created by Gugli on 14/08/2025.
//
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "integrators.h"
#include "raylib.h"

#define NORM2(x, y) ((x)*(x)+(y)*(y))

void DrawArrow(int startX, int startY, int endX, int endY, Color color) {
    float dx = endX - startX;
    float dy = endY - startY;
    float length = sqrtf(NORM2(dx, dy));

    // For the body, draw a line as thick as 8% of the body
    Vector2 bodyStart = {startX, startY};
    Vector2 bodyEnd = {endX - 0.2f*dx, endY - 0.2f*dy};
    float thickness = 0.05*length;
    DrawLineEx(bodyStart, bodyEnd, thickness, color);


    // Triangle vertices
    Vector2 p1 = {endX, endY}; // tip
    Vector2 p2 = {bodyEnd.x + 0.1f*dy, bodyEnd.y - 0.1f*dx};
    Vector2 p3 = {bodyEnd.x - 0.1f*dy, bodyEnd.y + 0.1f*dx};

    DrawTriangle(p1, p2, p3, color);
}


void harmonic_f(const PhysicsSystem* ps, double* a) {
    for (int i = 0; i < 3*ps->N; i+=3) {
        a[i] = -1*ps->x[i];
    }
}

int main(void) {
    const int screenWidth = 800;
    const int screenHeight = 450;

    printf("Init PhysicsSystem!\n");
    PhysicsSystem* ps = malloc(sizeof(PhysicsSystem));
    ps->N = 50;
    ps->m = calloc(3*ps->N, sizeof(*(ps->m)));
    ps->t = 0;
    ps->x = calloc(3*ps->N, sizeof(*(ps->x)));
    ps->v = calloc(3*ps->N, sizeof(*(ps->v)));

    for (int i = 0; i < ps->N; i++) {
        printf("Init particle %d\n", i);
        ps->m[i] = 1;
        ps->x[3*i] = cos( 2*M_PI*i/ps->N);
        ps->v[3*i] = sin( 2*M_PI*i/ps->N);
        ps->x[3*i+1] = (double) ( (screenWidth*(i+1))/(ps->N+1)); // Anchor point
        //printf("%lf\n", ps->x[3*i+1]);
        // Ignore z coordinate
    }

    //printf("Init Integrator!\n");
    Integrator* integrator = malloc(sizeof(Integrator));
    integrator->dt = 0.001;
    integrator->f = harmonic_f;
    init_RK4(integrator, ps);

    /*
    printf("Opening File\n");
    FILE* fp = fopen("output.csv", "w");

    printf("Starting simulation\n");
    for (int step = 0; step < 1e3; step++) {
         integrator->integrate(ps, integrator);
         if(step % 5 == 0) {
             fprintf(fp, "%f ", ps->x[0]);
             for(int i = 3; i < 3*ps->N; i+=3) {
                 fprintf(fp, ", %f", ps->x[i]);
             }
             fprintf(fp, "\n");
         }
    }
    fclose(fp);*/

    // Instead of saving the output to file, let's show it on a window with raylib
    InitWindow(screenWidth, screenHeight, "Observe this harmonic motion and relax");
    SetTargetFPS(60); // Set our game to run at 60 frames-per-second

    double energy;
    double ave_energy;
    double std_energy;
    char text[100];

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        for (int step=0; step<10; step++)
            integrator->integrate(ps, integrator);
        ave_energy = 0;
        std_energy = 0;

        for (int i = 0; i < ps->N; i++) {
            energy = 0.5*ps->v[3*i]*ps->v[3*i] + ps->x[3*i]*ps->x[3*i];
            ave_energy += energy;
            std_energy += (energy*energy); // Accumulate squares
            //printf("Energy of particle %d: %lf\n", i, energy);
        }
        ave_energy /= (double) ps->N;
        std_energy = sqrt(( std_energy/(double) ps->N - ave_energy*ave_energy) );

        //----------------------------------------------------------------------------------


        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
        ClearBackground(DARKBLUE);

        for (int i = 0; i < ps->N; i++) {
            int centerY = screenHeight/2 + (int) (0.50*screenHeight/2 *ps->x[3*i]);
            int centerX = (int) (ps->x[3*i+1]);
            int speedY = (int) (0.50*screenHeight/2 * ps->v[3*i]);
            float radius = 0.5*screenHeight/ps->N;

            DrawLine(centerX, screenHeight/2, centerX, centerY, BLACK);
            DrawCircle(centerX, centerY, radius, SKYBLUE);
            DrawArrow(centerX, centerY, centerX, centerY+speedY, SKYBLUE);
            DrawCircleLines(centerX, centerY, radius, BLACK);
        }

        sprintf(text, "Energy: %.3lf±%.3lf", ave_energy, std_energy);
        DrawText(text, 10, 50, 30, RAYWHITE);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // ----------------------- De-Initialization -------------------------------------------
    // Here I should dealloc eveything...
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------
    free_RK4(integrator);
    free(integrator);
    free(ps);


    return 0;
}