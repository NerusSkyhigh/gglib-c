#include <stdio.h>
#include <stdlib.h>

#include "msd.h"
//#include "rarray.h"
#include "parser.h"
#include "analysis.h"

#include <math.h>

int main(int argc, char* argv[]) {
    // Let's read filename and output_file from the command line
    if (argc != 5) {
        printf("Usage: %s <filename> <output_file_g1> <output_file_g2> <output_file_g3>\n", argv[0]);
        printf("Input provided:\n");
        for(int a=0; a<argc; a++) {
            printf("\t%s\n", argv[a]);
        }
        return 1;
    }
    const char* inp_file = argv[1];
    const char* output_file_g1 = argv[2];
    const char* output_file_g2 = argv[3];
    const char* output_file_g3 = argv[4];

    const int64_t T_EQ = 100000000;
    LAMMPSData data;

    printf("Reading LAMMPS data from %s\n", inp_file);
    printf("Setting T_EQ = %ld\n", T_EQ);
    initLAMMPSData(inp_file, &data, T_EQ);
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
    readLAMMPSCoordinates(inp_file, &data, T_EQ);
    printf("\tFirst coordinate: %f %f %f\n", data.coordinates[0], data.coordinates[1], data.coordinates[2]);

    int64_t offset = (data.num_timesteps-1)*(data.num_atoms)*3 + (data.num_atoms-1)*3;
    printf("\tLast coordinate:  %f %f %f\n", data.coordinates[offset], data.coordinates[offset+1], data.coordinates[offset+2]);


    //printf("Resaving file to output.lammpstrj as a tests\n");
    //writeLAMMPSData("output.lammpstrj", &data);


    /* -------------------------------------------------------------- */
    // Let's extract the relevant beads based on the molId
    //int64_t MOLID_MAXIRINGS[3] = { 605, 606, 607 };
    // TEMPLATE - (sim:LB_rev1)
    // Use notebook for minirings and number; find and replace spaces with `, ` to get the array
    //const int64_t NUMBER_OF_MOLECULES = 80;
    //int64_t MOLID_RELEVANT_MINIRINGS[80] = {7, 8, 9, 11, 16, 17, 18, 23, 24, 27, 30, 32, 40, 42, 54, 58, 59, 61, 62, 63,
    //                                       64, 80, 81, 86, 88, 89, 109, 111, 113, 142, 172, 173, 174, 175, 176, 178, 181, 183, 209, 216,
    //                                      217, 219, 246, 247, 252, 254, 282, 288, 289, 359, 362, 389, 395, 397, 422, 429, 430, 431, 457, 461,
    //                                      463, 486, 496, 497, 525, 528, 531, 532, 545, 553, 559, 562, 564, 565, 566, 577, 579, 581, 583, 585};

    // TEMPLATE - (sim:LB_rev0)
    //const int64_t NUMBER_OF_MOLECULES = 52;
    //int64_t MOLID_RELEVANT_MINIRINGS[52] = {7, 8, 23, 30, 50, 53, 54, 62, 64, 90, 113, 117, 142, 143, 144, 147, 148, 149, 174, 176,
    //                                      177, 181, 209, 210, 211, 212, 216, 217, 245, 246, 282, 288, 318, 319, 355, 356, 390, 431, 432, 456,
    //                                      495, 543, 552, 559, 566, 574, 575, 577, 578, 582, 585, 586};

    // Maxirings for all sims
    //const int64_t NUMBER_OF_MOLECULES = 3;
    //int64_t MOLID_RELEVANT_MINIRINGS[3] = { 605, 606, 607 };

    const int64_t NUMBER_OF_MOLECULES = 1;
    int64_t MOLID_RELEVANT_MINIRINGS[3] = { 7 };

    int64_t numOfSelectedAtoms;
    int64_t* sel_atIds = getAtomsFromMoleculeList(&data, MOLID_RELEVANT_MINIRINGS, NUMBER_OF_MOLECULES, &numOfSelectedAtoms);
    if (numOfSelectedAtoms < 1) {
        printf("No atom found. Closing.\n");
        exit(0);
    }
    printf("I found %ld atoms in the %ld mol specified. I expect %ld=num mol*60 atoms\n", numOfSelectedAtoms, NUMBER_OF_MOLECULES, NUMBER_OF_MOLECULES*60);

    /* -------------------------------------------------------------- */
    //const int64_t TARGET_ATOM_ID = 2;
    //printf("Extracting the atoms with type %ld.\n", TARGET_ATOM_ID);
    //int64_t numOfSelectedAtoms;
    //int64_t* sel_atIds = getAtomsFromType(&data, TARGET_ATOM_ID, &numOfSelectedAtoms);
    //if (numOfSelectedAtoms < 1) {
    //    printf("No atom found. Closing.\n");
    //    exit(0);
    //}
    //printf("I found %ld atoms with type %ld\n", numOfSelectedAtoms, TARGET_ATOM_ID);


    /* -------------------------------------------------------------- */
    printf("Computing g1 (Average MSD of all the beads)\n");

    // The maximum length possible is the last timestep divided by deltaTimestep rounded up
    int64_t total_n_windows = (data.timesteps[data.num_timesteps-1] - data.timesteps[0] + data.deltaTimestep - 1 ) / data.deltaTimestep;
    total_n_windows+=1; // Zero length window. I expect it to be zero and compute that as a check

    double* ave_g1 = (double*) calloc(total_n_windows, sizeof(double));

    // Average over particles
    for(int64_t p=0; p<numOfSelectedAtoms; p++) {
        // [TODO] Here I am shooting myself in the jewellery. If the code is too slow, this is one point to improve.
        // [TODO] Here I am shooting myself in the jewellery. If the code is too slow, this is one point to improve.
        double* bead_coord = getBeadTrajectory(&data, sel_atIds[p]);

        double* bead_msd   = compute_time_averaged_msd(bead_coord, data.timesteps, data.num_timesteps, data.deltaTimestep);

        // Accumulate the MSD into the average vector
        for(int64_t t=0; t<total_n_windows; t++) {
            ave_g1[t] += bead_msd[t] / ( (double) numOfSelectedAtoms);
        }
        free(bead_coord);
        free(bead_msd);
    }

    // Let's now print the output to file
    FILE* output_file = fopen(output_file_g1, "w");
    fprintf(output_file, "# TIMESTEP ave_g1\n");
    for (int64_t t=0; t<total_n_windows; t++) {
        fprintf(output_file, "%10ld %12.6f\n", t*data.deltaTimestep, ave_g1[t]);
    }
    fclose(output_file);
    free(ave_g1);


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
    }


    /* -------------------------------------------------------------- */
    printf("Removing center of mass motion\n");
    for(int64_t t=0; t<data.num_timesteps; t++) {
        int64_t offset = data.num_atoms*t*3;
        for(int64_t i=0; i<data.num_atoms; i++) {
            data.coordinates[offset+3*i]   -= com[t*3];
            data.coordinates[offset+3*i+1] -= com[t*3+1];
            data.coordinates[offset+3*i+2] -= com[t*3+2];
            //printf("%lf %lf %lf\t", com[t*3], com[t*3+1], com[t*3+2]);
            //printf("%lf %lf %lf\n", data.coordinates[t*3], data.coordinates[t*3+1], data.coordinates[t*3+2]);
        }
    }

    printf("Computing g2 (Average MSD of all the beads)\n");
    double* ave_g2 = (double*) calloc(total_n_windows, sizeof(double));
    for(int64_t p=0; p<numOfSelectedAtoms; p++) {
        double* bead_coord = getBeadTrajectory(&data, sel_atIds[p]);
        double* bead_msd   = compute_time_averaged_msd(bead_coord, data.timesteps, data.num_timesteps, data.deltaTimestep);
        for(int64_t t=0; t<total_n_windows; t++) {
            ave_g2[t] += bead_msd[t]/ ( (double) numOfSelectedAtoms);
        }
        free(bead_coord);
        free(bead_msd);
    }

    // Let's now print the output to file
    output_file = fopen(output_file_g2, "w");
    fprintf(output_file, "# TIMESTEP ave_g2\n");
    for (int64_t t=0; t<total_n_windows; t++) {
        fprintf(output_file, "%10ld %12.6f\n", t*data.deltaTimestep, ave_g2[t]);
    }
    fclose(output_file);
    free(ave_g2);

    /* -------------------------------------------------------------- */
    printf("Computing g3 (Average MSD of the com)\n");
    double* ave_g3 = (double*) calloc(total_n_windows, sizeof(double));
    double* com_msd   = compute_time_averaged_msd(com, data.timesteps, data.num_timesteps, data.deltaTimestep);
    for(int64_t t=0; t<total_n_windows; t++) {
        ave_g3[t] += com_msd[t];
    }
    free(com_msd);

    // Let's now print the output to file
    output_file = fopen(output_file_g3, "w");
    fprintf(output_file, "# TIMESTEP ave_g3\n");
    for (int64_t t=0; t<total_n_windows; t++) {
        fprintf(output_file, "%10ld %12.6f\n", t*data.deltaTimestep, ave_g3[t]);
    }
    fclose(output_file);
    free(ave_g3);




    freeLAMMPSData(&data);

    free(sel_atIds);
    free(com);

    return 0;
}
