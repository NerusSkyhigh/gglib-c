#ifndef LAMMPS_DATA_H
#define LAMMPS_DATA_H

#include <stdint.h>

// Struct to hold LAMMPS data
typedef struct {
    // TODO: Have the struct contain a copy of the topology
    //       so that it becomes an independent object.
    int64_t num_atoms;        // Number of atoms
    double* box;              // Box dimensions
    int64_t* atomIds;         // Atom IDs array
    int64_t* moleculeIds;     // Molecule IDs array
    int64_t* atomTypes;       // Atom types array
    int64_t timestep;         // Timestep of the frame
    int64_t index;            // Index of the timestep in the

    double* coordinates;      // Coordinates array (x, y, z for each atom)
} LAMMPSFrame;


// Struct to hold LAMMPS data
typedef struct {
    int64_t deltaTimestep;    // Delta timestep
    int64_t num_atoms;        // Number of atoms
    int64_t num_timesteps;    // Number of timesteps
    double* box;              // Box dimensions
    int64_t* atomIds;         // Atom IDs array
    int64_t* moleculeIds;     // Molecule IDs array
    int64_t* atomTypes;       // Atom types array
    //double* coordinates;      // Coordinates array (x, y, z for each atom)
    int64_t* timesteps;       // Timesteps array
    long* frame_offsets;      // Record the file offset of each frame for lazy load
                              // Note that it points to the line after ITEM: TIMESTEP
    LAMMPSFrame** frames;
} LAMMPSData;



// Function declarations
void initLAMMPSData(const char* filename, LAMMPSData* data, const int64_t T_EQ);
void read_frame_at_offset(FILE* file, LAMMPSFrame* frame, long offset);
void readLAMMPSCoordinates(const char* filename, LAMMPSData* data, const int64_t T_EQ);

// Helper functions to write TIMESTEPS, BOX, N_ATOMS and ATOMS information
void writeLAMMPSData(const char* filename, const LAMMPSData* data);

// Free functions
void freeLAMMPSFrame(LAMMPSFrame* frame);
void freeLAMMPSData(LAMMPSData* data);

// Check functions
void checkTimestepMismatch(const LAMMPSData* data);

#endif // LAMMPS_DATA_H
