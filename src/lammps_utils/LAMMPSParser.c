#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rarray.h"
#include "LAMMPSParser.h"

#include "fileutils.h"
#include "memutils.h"

//**************** READER FUNCTION ****************
// WARNING: The reader functions assume that the
// file pointer is located after the corresponding
// header file the writer functions assume that the
// file pointer is located before the header
void initLAMMPSData(const char* filename, LAMMPSData* data, const int64_t T_EQ) {
    FILE* file;
    SAFE_FILE_OPEN(file, filename, "r", "");

    char line[4096];
    int64_t timestep,
            atomId, molId, atomType;
    int64_t natoms_found = 0;
    int64_t box_found = 0;
    int64_t atom_list_constructed = 0;
    double box_min, box_max;

    // Temporary buffers are constructed using rarray
    size_t initialCapacity = 1000;
    rarray* timesteps_buf   = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* atomIds_buf     = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* moleculeIds_buf = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* atomTypes_buf   = rarray_init(sizeof(int64_t), initialCapacity);
    rarray* box_buff        = rarray_init(sizeof(double),  initialCapacity);
    rarray* fo_buff         = rarray_init(sizeof(long), initialCapacity);

    while(fgets(line, sizeof(line), file)) {
        // The first time we see the atom list, we record it.
        // NOTE: if ITEM strings are reused across the codebase, consider putting them in a #define
        //#define ITEM_HEADER "ITEM: ATOMS id mol type xu yu zu"
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
            long offset = ftell(file);
            if (fgets(line, sizeof(line), file)) {
                timestep = atoll(line);

                // Add the timestep to the list only if it is beyond the equilibration time
                if(timestep >= T_EQ) {
                    rarray_push(timesteps_buf, &timestep);
                    rarray_push(fo_buff, &offset);
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
    SAFE_FILE_CLOSE(file);
    //fclose(file);

    // Convert rarray buffers to normal arrays and store them in the LAMMPSData struct
    {
        data->num_timesteps = (int64_t)     rarray_size(timesteps_buf);
        data->timesteps     = (int64_t*)    rarray_to_array(timesteps_buf);
        data->atomIds       = (int64_t*)    rarray_to_array(atomIds_buf);
        data->moleculeIds   = (int64_t*)    rarray_to_array(moleculeIds_buf);
        data->atomTypes     = (int64_t*)    rarray_to_array(atomTypes_buf);
        data->box           = (double*)     rarray_to_array(box_buff);
        data->frame_offsets = (long*)       rarray_to_array(fo_buff);
    }

    if (rarray_size(timesteps_buf) != rarray_size(fo_buff)) {
        fprintf(stderr, "Internal error: timestep/offset mismatch\n");
        exit(1);
    }

    // Free rarray structures (the data buffers are now owned by LAMMPSData)
    {
        rarray_free(timesteps_buf);
        rarray_free(atomIds_buf);
        rarray_free(moleculeIds_buf);
        rarray_free(atomTypes_buf);
        rarray_free(box_buff);
        rarray_free(fo_buff);
    }

    // Compute deltaTimestep
    if (data->num_timesteps > 1) {
        data->deltaTimestep = data->timesteps[1] - data->timesteps[0];
    } else {
        fprintf(stderr, "Error: less than 2 timesteps found.\n");
        data->deltaTimestep = 0;
    }

}

void read_frame_at_offset(FILE* file, LAMMPSFrame* frame, const long offset) {
    char line[1024];
    int64_t atomId, molId, atomType;
    double x, y, z;

    // Move the file pointer to the correct position and read the timestep value
    fseek(file, offset, SEEK_SET);

    // Read the timestep and check it is the correct one
    if (fgets(line, sizeof(line), file)) {
        int64_t timestep = atoll(line);
        if (timestep != frame->timestep) {
            fprintf(stderr, "Internal error: timestep/offset-timestep mismatch [%ld/%ld]\n",
                        frame->timestep, timestep);
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "Error reading frame at %ld\n", offset);
        exit(EXIT_FAILURE);
    }

    // Read data to reach coordinate section
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "ITEM: ATOMS id mol type xu yu zu")) {
            break;
        }
    }

    // Store data into the frame
    char errorMsg[4096];
    snprintf(errorMsg,4096, "Unable to allocate memory for frame [timestep %ld; offset=%ld]", frame->timestep, offset);
    SAFE_MALLOC(frame->coordinates, 3 * frame->num_atoms * sizeof(double), errorMsg);
    // OLD: frame->coordinates = (double*) malloc(3 * frame->num_atoms * sizeof(double));

    int64_t a_idx = 0;
    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "%ld %ld %ld %lf %lf %lf", &atomId, &molId, &atomType, &x, &y, &z) == 6) {
            frame->coordinates[a_idx] = x;
            frame->coordinates[a_idx + 1] = y;
            frame->coordinates[a_idx + 2] = z;
            a_idx +=3;
            //TODO: Add support for unsorted arrays. Right now this check
            //      is valid only if the atoms are sorted
            if (atomId == frame->atomIds[frame->num_atoms-1]) {
                break;
            }
        }
    }
    if (a_idx != 3*frame->num_atoms) {
        fprintf(stderr, "Internal error: expected %ld atoms but found %ld\n", frame->num_atoms, a_idx);
        exit(EXIT_FAILURE);
    }

}

void readLAMMPSCoordinates(const char* filename, LAMMPSData* data, const int64_t T_EQ) {
    FILE* file;
    SAFE_FILE_OPEN(file, filename, "r", "");


    char line[1024];
    int atomId, molId, atomType;
    double x, y, z;

    // Let's allocate a vector for the frames
    SAFE_MALLOC(data->frames, data->num_timesteps* sizeof(LAMMPSFrame*), "Unable to allocate memory for frames array");

    int64_t timestep = 0;

    for (int64_t frameIdx=0; frameIdx < data->num_timesteps; frameIdx++) {
        char errorMsg[4096];
        snprintf(errorMsg,4096, "Unable to allocate memory for frame struct [timestep %ld; offset=%ld]", data->timesteps[frameIdx], data->frame_offsets[frameIdx]);
        SAFE_MALLOC(data->frames[frameIdx], sizeof(LAMMPSFrame), errorMsg);
        data->frames[frameIdx]->timestep = data->timesteps[frameIdx];
        data->frames[frameIdx]->num_atoms = data->num_atoms;
        data->frames[frameIdx]->atomIds = data->atomIds;
        data->frames[frameIdx]->atomTypes = data->atomTypes;
        data->frames[frameIdx]->moleculeIds = data->moleculeIds;
        data->frames[frameIdx]->box = data->box;

        read_frame_at_offset(file, data->frames[frameIdx], data->frame_offsets[frameIdx]);
    }
    fclose(file);

}

//**************** WRITER FUNCTIONS ****************
// WARNING: The reader functions assume that the
// file pointer is located after the corresponding
// header file the writer functions assume that the
// file pointer is located before the header
static inline void writeTimestep(FILE* file, LAMMPSFrame* frame) {
    fprintf(file, "ITEM: TIMESTEP\n");
    fprintf(file, "%ld\n", frame->timestep);
}
static inline void writeNAtoms(FILE* file, LAMMPSFrame* frame) {
    fprintf(file, "ITEM: NUMBER OF ATOMS\n");
    fprintf(file, "%ld\n", frame->num_atoms);
}
static inline void writeBox(FILE* file, LAMMPSFrame* frame) {
    fprintf(file, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(file, "%.16e %.16e\n", frame->box[0], frame->box[1]);
    fprintf(file, "%.16e %.16e\n", frame->box[2], frame->box[3]);
    fprintf(file, "%.16e %.16e\n", frame->box[4], frame->box[5]);
}

static inline void writeCoordinates(FILE* file, LAMMPSFrame* frame) {
    fprintf(file, "ITEM: ATOMS id mol type xu yu zu\n");
    for (int i = 0; i < frame->num_atoms; i++) {
        fprintf(file, "%10ld %6ld %3ld %12.6f %12.6f %12.6f\n",
                frame->atomIds[i], frame->moleculeIds[i], frame->atomTypes[i],
                frame->coordinates[3*i], frame->coordinates[3*i+1], frame->coordinates[3*i+2]);
    }
}


void writeLAMMPSData(const char* filename, const LAMMPSData* data) {
    FILE* file;
    SAFE_FILE_OPEN(file, filename, "w", "");

    // Write timesteps one at a time
    for (int64_t t = 0; t < data->num_timesteps; t++) {
        LAMMPSFrame* current_frame = data->frames[t];
        writeTimestep(file, current_frame);
        writeNAtoms(file, current_frame);
        writeBox(file, current_frame);
        writeCoordinates(file, current_frame);
    }

    fclose(file);
}


//**************** FREE STRUCT FUNCTIONS ****************
void freeLAMMPSFrame(LAMMPSFrame* frame) {
    SAFE_FREE(frame->coordinates);
}

void freeLAMMPSData(LAMMPSData* data) {
    SAFE_FREE(data->timesteps);
    SAFE_FREE(data->atomIds);
    SAFE_FREE(data->moleculeIds);
    SAFE_FREE(data->atomTypes);
    SAFE_FREE(data->box);
    for (int64_t frameIdx=0; frameIdx < data->num_timesteps; frameIdx++) {
        freeLAMMPSFrame(data->frames[frameIdx]);
        SAFE_FREE(data->frames[frameIdx]);
    }
    SAFE_FREE(data->frames);

}

//**************** CHECK FUNCTIONS ****************
void checkTimestepMismatch(const LAMMPSData* data) {
    printf("Checking for timestep mismatch\n");
    int64_t num_timesteps = data->num_timesteps;

    for (int64_t i = 1; i < num_timesteps; i++) {
        int64_t current  = data->timesteps[i];
        int64_t previous = data->timesteps[i-1];

        if( (current - previous) != data->deltaTimestep) {
            fprintf(stderr, "Timestep mismatch at index %ld\n", i);
            fprintf(stderr, "Previous timestep: %ld\n", previous);
            fprintf(stderr, "Current timestep: %ld\n", current);
            fprintf(stderr, "Delta timestep: %ld\n", data->deltaTimestep);
            fprintf(stderr, "Delta index: %ld\n", (current - previous) / data->deltaTimestep);
        }
    }
}