//
// Created by gu on 26/03/25.
//
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>

#include "LAMMPSParser.h"
#include "LAMMPSDataAnalysis.h"

int main(int argc, char* argv[]) {
    // Let's read filename and output_file from the command line
    if (argc != 4) {
        printf("Usage: %s <filename> <output_file_com> <output_file_com_removed>\n", argv[0]);
        return 1;
    }
    const char* inp_file = argv[1];
    const char* output_file_com = argv[2];
    const char* output_file_traj_com_removed = argv[3];

    LAMMPSData data;

    printf("Reading LAMMPS data from %s\n", inp_file);
    initLAMMPSData(inp_file, &data);
    checkTimestepMismatch(&data);

    printf("\tNumber of atoms: %ld\n", data.num_atoms);
    printf("\tNumber of timesteps: %ld\n", data.num_timesteps);
    printf("\tDelta timestep: %ld\n", data.deltaTimestep);
    printf("\tFirst timestep: %ld\n", data.timesteps[0]);
    printf("\tLast timestep: %ld\n", data.timesteps[data.num_timesteps-1]);
    printf("\tFirst atom: %ld\n", data.atomIds[0]);
    printf("\tFirst molecule: %ld\n", data.moleculeIds[0]);
    printf("\tFirst atom type: %ld\n", data.atomTypes[0]);

    printf("BOX:\n");
    printf("\t%.16e %.16e\n", data.box[0], data.box[1]);
    printf("\t%.16e %.16e\n", data.box[2], data.box[3]);
    printf("\t%.16e %.16e\n", data.box[4], data.box[5]);

    printf("COORDINATES:\n");
    readLAMMPSCoordinates(inp_file, &data);
    printf("\tFirst coordinate: %f %f %f\n", data.coordinates[0], data.coordinates[1], data.coordinates[2]);
    printf("\tLast coordinate:  %f %f %f\n", data.coordinates[3*data.num_atoms-3], data.coordinates[3*data.num_atoms-2], data.coordinates[3*data.num_atoms-1]);


    //printf("Resaving file to output.lammpstrj as a tests\n");
    //writeLAMMPSData("output.lammpstrj", &data);

    /* -------------------------------------------------------------- */
    printf("Computing CoM coordinates\n");
    double* com = (double*) malloc(sizeof(double)*data.num_timesteps*3);

    for(int64_t t=0; t<data.num_timesteps; t++) {
        // As the com array has structure [x0 y0 z0 x1 y1 z1 ... xn yn zn],
        // `com+t*3` points to xt.
        // Similarly, `data.coordinates + data.num_atoms*t*3` points to
        // the x coordinate of the first particle at time t.
        size_t offset = (size_t) data.num_atoms*t*3;
        compute_CoM(data.coordinates+offset, data.num_atoms, com+3*t);

        //fprintf(stderr, "%10ld %10ld %12.6lf %12.6lf %12.6lf [%ld]\n", t, data.timesteps[t], com[0], com[1], com[2], offset);
    }

    /* -------------------------------------------------------------- */
    printf("Removing center of mass motion\n");
    for(int64_t t=0; t<data.num_timesteps; t++) {
        int64_t offset = data.num_atoms*t*3;
        for(int64_t i=0; i<data.num_atoms; i++) {
            data.coordinates[offset+3*i]   -= com[t*3];
            data.coordinates[offset+3*i+1] -= com[t*3+1];
            data.coordinates[offset+3*i+2] -= com[t*3+2];
        }
    }

    /* -------------------------------------------------------------- */
    // Let's also save the com to a file
    printf("Saving com data to %s\n", output_file_com);
    FILE* file = fopen(output_file_com, "w");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", output_file_com);
    } else {
        for(int64_t t=0; t<data.num_timesteps; t++) {
            fprintf(file, "%10ld %12.6lf %12.6lf %12.6lf\n", data.timesteps[t], com[t*3], com[t*3+1], com[t*3+2]);
        }
        fclose(file);
    }

    writeLAMMPSData(output_file_traj_com_removed, &data);

    freeLAMMPSData(&data);
    free(com);

    return 0;
}
