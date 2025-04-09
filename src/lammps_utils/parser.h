#ifndef LAMMPS_DATA_H
#define LAMMPS_DATA_H

#include <stdint.h>

// Struct to hold LAMMPS data
typedef struct {
    int64_t deltaTimestep;    // Delta timestep
    int64_t num_atoms;        // Number of atoms
    int64_t num_timesteps;    // Number of timesteps
    double* box;              // Box dimensions
    int64_t* atomIds;         // Atom IDs array
    int64_t* moleculeIds;     // Molecule IDs array
    int64_t* atomTypes;       // Atom types array
    double* coordinates;      // Coordinates array (x, y, z for each atom)
    int64_t* timesteps;       // Timesteps array
} LAMMPSData;

// Function declarations
void initLAMMPSData(const char* filename, LAMMPSData* data, const int64_t T_EQ);
void readLAMMPSCoordinates(const char* filename, LAMMPSData* data, const int64_t T_EQ);
void checkTimestepMismatch(const LAMMPSData* data);
void freeLAMMPSData(LAMMPSData* data);
void writeLAMMPSData(const char* filename, const LAMMPSData* data);

double* getBeadTrajectory(const LAMMPSData* data, const int64_t atom_id);
int64_t* getAtomsFromMoleculeList(const LAMMPSData* data, const int64_t* molecule_ids, const int64_t numMolecules, int64_t* numOfSelectedAtoms);
int64_t* getAtomsInMolecule(const LAMMPSData* data, const int64_t molecule_id, int64_t* numOfSelectedAtoms);
int64_t* getAtomsFromType(const LAMMPSData* data, const int64_t a_type, int64_t* numOfSelectedAtoms);
#endif // LAMMPS_DATA_H
