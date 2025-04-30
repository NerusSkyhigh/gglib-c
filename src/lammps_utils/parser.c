#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rarray.h"
#include "parser.h"

void initLAMMPSData(const char* filename, LAMMPSData* data, const int64_t T_EQ) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    char line[4096];
    int64_t timestep;
    int64_t atomId, molId, atomType;
    int64_t natoms_found = 0;
    int box_found = 0;
    int atom_list_constructed = 0;
    double box_min, box_max;


    // Temporary buffers using rarray
    size_t initialCapacity = 1; //1000;
    rarray* timesteps_buf   = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* atomIds_buf     = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* moleculeIds_buf = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* atomTypes_buf   = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* box_buff        = rarray_init(sizeof(double),  initialCapacity);

    while(fgets(line, sizeof(line), file)) {
        // The first time we see the atom list, we construct the atom list
        if (!atom_list_constructed && strstr(line, "ITEM: ATOMS id mol type xu yu zu")) {
            while(fgets(line, sizeof(line), file)) {
                if (strstr(line, "ITEM")) {
                    break;
                }
                if (sscanf(line, "%ld %ld %ld", &atomId, &molId, &atomType) == 3) {
                    rarray_push(atomIds_buf, &atomId);
                    rarray_push(moleculeIds_buf, &molId);
                    rarray_push(atomTypes_buf, &atomType);
                }
            }
            atom_list_constructed = 1;
        }

        // Let's also read the box
        if (!box_found && strstr(line, "ITEM: BOX BOUNDS pp pp pp")) {
            while(fgets(line, sizeof(line), file)) {
                if (strstr(line, "ITEM")) {
                    break;
                }
                if (sscanf(line, "%lf %lf", &box_min, &box_max) == 2) {
                    rarray_push(box_buff, &box_min);
                    rarray_push(box_buff, &box_max);
                }
            }
            box_found = 1;
        }


        // Timesteps are located after the `ITEM: TIMESTEP` line
        if (strstr(line, "ITEM: TIMESTEP")) {
            if (fgets(line, sizeof(line), file)) {
                timestep = atoll(line);

                // Add the timestep to the list only if it is beyond the equilibration time
                if(timestep >= T_EQ) {
                    rarray_push(timesteps_buf, &timestep);
                }
            }
        }
        
        if (!natoms_found && strstr(line, "ITEM: NUMBER OF ATOMS")) {
            if (fgets(line, sizeof(line), file)) {
                data->num_atoms = atoi(line);
                natoms_found = 1;
            }
        }
    }
    fclose(file);

    // Convert rarray buffers to normal arrays
    data->num_timesteps = (int64_t)     rarray_size(timesteps_buf);
    data->timesteps     = (int64_t*)    rarray_to_array(timesteps_buf);
    data->atomIds       = (int64_t*)    rarray_to_array(atomIds_buf);
    data->moleculeIds   = (int64_t*)    rarray_to_array(moleculeIds_buf);
    data->atomTypes     = (int64_t*)    rarray_to_array(atomTypes_buf);
    data->box           = (double*)     rarray_to_array(box_buff);

    // Free rarray structures (the data buffers are now owned by LAMMPSData)
    rarray_free(timesteps_buf);
    rarray_free(atomIds_buf);
    rarray_free(moleculeIds_buf);
    rarray_free(atomTypes_buf);
    rarray_free(box_buff);

    // Compute deltaTimestep
    if (data->num_timesteps > 1) {
        data->deltaTimestep = data->timesteps[1] - data->timesteps[0];
    } else {
        fprintf(stderr, "Error: less than 2 timesteps found.\n");
        data->deltaTimestep = 0;
    }
}



void readLAMMPSCoordinates(const char* filename, LAMMPSData* data, const int64_t T_EQ) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    char line[1024];
    int atomId, molId, atomType;
    double x, y, z;

    rarray* coordinates_buff   = rarray_init(sizeof(double), 10);

    int64_t timestep = 0;

    while (fgets(line, sizeof(line), file)) {

        if (strstr(line, "ITEM: TIMESTEP")) {
            if (fgets(line, sizeof(line), file)) {
                timestep = atoll(line);
            }
        }

        if (strstr(line, "ITEM: ATOMS id mol type xu yu zu")) {
            while (fgets(line, sizeof(line), file)) {
                if (sscanf(line, "%d %d %d %lf %lf %lf", &atomId, &molId, &atomType, &x, &y, &z) == 6) {
                    if (timestep >= T_EQ) {
                        rarray_push(coordinates_buff, &x);
                        rarray_push(coordinates_buff, &y);
                        rarray_push(coordinates_buff, &z);
                    }
                    if (atomId == data->atomIds[data->num_atoms-1]) {
                        break;
                    }
                }
            }
        }

    }
    fclose(file);

    data->coordinates = (double*) rarray_to_array(coordinates_buff);
    rarray_free(coordinates_buff);
}

void freeLAMMPSData(LAMMPSData* data) {
    free(data->timesteps);
    free(data->atomIds);
    free(data->coordinates);
    free(data->moleculeIds);
    free(data->atomTypes);
    free(data->box);
}

void checkTimestepMismatch(const LAMMPSData* data) {
    printf("Checking for timestep mismatch\n");
    int64_t num_timesteps = data->num_timesteps;

    for (int64_t i = 1; i < num_timesteps; i++) {
        int64_t current  = data->timesteps[i];
        int64_t previous = data->timesteps[i-1];

        if( (current - previous) != data->deltaTimestep) {
            fprintf(stderr, "Timestep mismatch at index %zu\n", i);
            fprintf(stderr, "Previous timestep: %ld\n", previous);
            fprintf(stderr, "Current timestep: %ld\n", current);
            fprintf(stderr, "Delta timestep: %ld\n", data->deltaTimestep);
            fprintf(stderr, "Delta index: %ld\n", (current - previous) / data->deltaTimestep);
        }
    }
}


void writeLAMMPSData(const char* filename, const LAMMPSData* data) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    // Write a single timestep
    for (int64_t t = 0; t < data->num_timesteps; t++)
    {
        int64_t offset = 3*t * data->num_atoms;
        fprintf(file, "ITEM: TIMESTEP\n");
        fprintf(file, "%ld\n", data->timesteps[t]);

        fprintf(file, "ITEM: NUMBER OF ATOMS\n");
        fprintf(file, "%ld\n", data->num_atoms);

        fprintf(file, "ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(file, "%.16e %.16e\n", data->box[0], data->box[1]);
        fprintf(file, "%.16e %.16e\n", data->box[2], data->box[3]);
        fprintf(file, "%.16e %.16e\n", data->box[4], data->box[5]);

        fprintf(file, "ITEM: ATOMS id mol type xu yu zu\n");
        for (int i = 0; i < data->num_atoms; i++) {
            fprintf(file, "%10ld %6ld %3ld %12.6f %12.6f %12.6f\n", data->atomIds[i], data->moleculeIds[i], data->atomTypes[i],
                    data->coordinates[offset+3*i], data->coordinates[offset+3*i+1], data->coordinates[offset+3*i+2]);
        }
    }

    fclose(file);
}