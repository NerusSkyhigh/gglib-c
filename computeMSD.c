#include <stdio.h>
#include <stdlib.h>

#include "msd.h"
#include "LAMMPSParser.h"
#include "analysis.h"
#include "configParser.h"
#include "rarray.h"


int main(const int argc, char* argv[]) {
    // Let's read filename and output_file from the command line
    if (argc != 2) {
        printf("Usage: %s <filename> <config.ini>\n", argv[0]);
        printf("Input provided:\n");
        for(int a=0; a<argc; a++) {
            printf("\t%s\n", argv[a]);
        }
        return 1;
    }
    const char* configfile = argv[1];

    // READ CONFIGURATION FILE
    struct config cfg;
    read_config(configfile, &cfg);

    printf("Config Loaded:\n");
    printf("\tT_EQ=%ld\n", cfg.T_EQ);

    printf("\tatom_ids=[ ");
    for(int64_t a=0; a<cfg.atom_ids_size; a++) {printf("%ld ", cfg.atom_ids[a]);}
    printf("]\n");
    printf("\tatom_ids_size=%ld\n", cfg.atom_ids_size);

    printf("\tmol_ids=[ ");
    for(int64_t m=0; m<cfg.mol_ids_size; m++) {printf("%ld ", cfg.mol_ids[m]);}
    printf("]\n");
    printf("\tmol_ids_size=%ld\n", cfg.mol_ids_size);

    printf("\tatom_ids=[ ");
    for(int64_t i=0; i<cfg.atom_types_size; i++) { printf("%ld ", cfg.atom_types[i]); }
    printf("]\n");
    printf("\tatom_ids_size=%ld\n", cfg.atom_ids_size);

    printf("\tinput_file=%s\n", cfg.input_file);
    printf("\tg1_output_file=%s\n", cfg.g1_output_file);
    printf("\tg2_output_file=%s\n", cfg.g2_output_file);
    printf("\tg3_output_file=%s\n", cfg.g3_output_file);
    printf("\tcom_output_file=%s\n", cfg.com_output_file);
    printf("\n");


    // READ LAMMPS DATA
    LAMMPSData data;

    printf("Reading LAMMPS data from %s\n", cfg.input_file);
    printf("Setting T_EQ = %ld\n", cfg.T_EQ);
    initLAMMPSData(cfg.input_file, &data, cfg.T_EQ);
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
    readLAMMPSCoordinates(cfg.input_file, &data, cfg.T_EQ);
    printf("\tFirst coordinate: %f %f %f\n", data.frames[0]->coordinates[0], data.frames[0]->coordinates[1], data.frames[0]->coordinates[2]);

    int64_t lastFrame = (data.num_timesteps-1);
    int64_t lastParticle = (data.num_atoms-1)*3;
    printf("\tLast coordinate:  %f %f %f\n", data.frames[lastFrame]->coordinates[lastParticle  ],
                                                 data.frames[lastFrame]->coordinates[lastParticle+1],
                                                 data.frames[lastFrame]->coordinates[lastParticle+2]);


    //printf("Resaving file to output.lammpstrj as a tests\n");
    //writeLAMMPSData("output.lammpstrj", &data);


    {;
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
    }

    // EXTRACT THE RELEVANT BEADS FROM THE TRAJECTORY
    size_t initialCapacity = 1;
    rarray* beads_buf   = rarray_init(sizeof(int64_t), initialCapacity);

    // Copy atom id of beads from the config file
    if (cfg.atom_ids_size > 0) {
        for(int64_t b=0; b<cfg.atom_ids_size; b++) {
            rarray_push(beads_buf, cfg.atom_ids+b);
        }
    } else {
        fprintf(stderr, "No bead list provided in the file.\n");
    }

    /* -------------------------------------------------------------- */
    // Select atom_id of beads with molId specified
    int64_t numOfSelectedAtoms;
    int64_t* sel_atIds;

    if (cfg.mol_ids_size > 0) {
        sel_atIds = getAtomsFromMoleculeList(&data, cfg.mol_ids, cfg.mol_ids_size, &numOfSelectedAtoms);
        if (numOfSelectedAtoms < 1) {
            printf("No atom found in molids specified in the file %s.\n", cfg.input_file);
        } else {
            for(int64_t b=0; b<numOfSelectedAtoms; b++) {
                rarray_push(beads_buf, sel_atIds+b);
            }
        }
        free(sel_atIds);
    }

    // Select atom_id of beads with type specified
    sel_atIds = getAtomsFromTypeList(&data, cfg.atom_types, cfg.atom_types_size, &numOfSelectedAtoms);
    if (numOfSelectedAtoms < 1) {
        printf("No atom of type specified in the file %s.\n", cfg.input_file);
    } else {
        // Copy atoms
        for(int64_t b=0; b<numOfSelectedAtoms; b++) {
            rarray_push(beads_buf, sel_atIds+b);
        }
    }
    free(sel_atIds);


    numOfSelectedAtoms =  rarray_size(beads_buf);
    sel_atIds = (int64_t*) rarray_to_array(beads_buf);
    rarray_free(beads_buf);
    if (numOfSelectedAtoms>0) {
        printf("I selected %ld of %ld atoms.\n", numOfSelectedAtoms, data.num_atoms);
    } else {
        fprintf(stderr, "No atom found in the file %s. Closing.\n", cfg.input_file);
        exit(-1);
    }


    /* -------------------------------------------------------------- */
    printf("Computing g1 (Average MSD of all the beads)\n");

    // The maximum length possible is the last timestep divided by deltaTimestep rounded up
    int64_t total_n_windows = (data.timesteps[data.num_timesteps-1] - data.timesteps[0] + data.deltaTimestep - 1 ) / data.deltaTimestep;
    total_n_windows+=1; // Zero length window. I expect it to be zero and compute that as a check

    double* ave_g1 = (double*) calloc(total_n_windows, sizeof(double));

    // Average over particles
    for(int64_t p=0; p<numOfSelectedAtoms; p++) {
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
    FILE* output_file = fopen(cfg.g1_output_file, "w");
    printf("\t saving to %s\n", cfg.g1_output_file);
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

        //OLD: size_t offset = (size_t) data.num_atoms*t*3;
        //OLD: compute_CoM(data.coordinates+offset, data.num_atoms, com+3*t);
        compute_CoM(data.frames[t]->coordinates, data.num_atoms, com+3*t);
    }

    output_file = fopen(cfg.com_output_file, "w");
    printf("\t saving to %s\n", cfg.com_output_file);
    fprintf(output_file, "# TIMESTEP x y z\n");
    for (int64_t t=0; t<data.num_timesteps; t++) {
        double x = com[3*t];
        double y = com[3*t+1];
        double z = com[3*t+2];
        fprintf(output_file, "%10ld %12.6f %12.6f %12.6f\n", t*data.deltaTimestep, x, y, z);
    }
    fclose(output_file);


    /* -------------------------------------------------------------- */
    printf("Removing center of mass motion\n");
    for(int64_t t=0; t<data.num_timesteps; t++) {
        //int64_t offset = data.num_atoms*t*3;
        for(int64_t i=0; i<data.num_atoms; i++) {
            data.frames[t]->coordinates[3*i]   -= com[t*3];     // OLD: data.coordinates[offset+3*i]   -= com[t*3];
            data.frames[t]->coordinates[3*i+1] -= com[t*3+1];   // OLD: data.coordinates[offset+3*i+1] -= com[t*3+1];
            data.frames[t]->coordinates[3*i+2] -= com[t*3+2];   // OLD: data.coordinates[offset+3*i+2] -= com[t*3+2];
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
    output_file = fopen(cfg.g2_output_file, "w");
    printf("\t saving to %s\n", cfg.g2_output_file);
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
    output_file = fopen(cfg.g3_output_file, "w");
    printf("\t saving to %s\n", cfg.g3_output_file);
    fprintf(output_file, "# TIMESTEP ave_g3\n");
    for (int64_t t=0; t<total_n_windows; t++) {
        fprintf(output_file, "%10ld %12.6f\n", t*data.deltaTimestep, ave_g3[t]);
    }
    fclose(output_file);
    free(ave_g3);




    freeLAMMPSData(&data);

    free(sel_atIds);
    free(com);
    free_config(&cfg);

    return 0;
}
