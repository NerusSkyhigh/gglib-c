//
// Created by gu on 27/03/25.
//
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "analysis.h"
#include "LAMMPSParser.h"
#include "rarray.h"

void compute_CoM(const double* frame, const int64_t num_atoms, double* com_coord) {
    // Set to zero for accumulation
    com_coord[0] = 0.;
    com_coord[1] = 0.;
    com_coord[2] = 0.;

    for (int i = 0; i < num_atoms; i++) {
        // Let's divide by number of atoms to avoid overflow
        com_coord[0] += frame[3l*i]   / ( (double) num_atoms);
        com_coord[1] += frame[3l*i+1l] / ( (double) num_atoms);
        com_coord[2] += frame[3l*i+2l] / ( (double) num_atoms);
    }
    //printf("%lf %lf %lf\n", com_coord[0], com_coord[1], com_coord[2]);
}

double* getBeadTrajectory(const LAMMPSData* data, const int64_t atom_id) {
    /*
     * Return the trajectory of a select bead as a function of time.
     * The array is a linear array.
     */
    rarray* coordinates_buff   = rarray_init(sizeof(double), 10);

    for (int64_t t = 0; t < data->num_timesteps; t++) {
        // Frame offset (elapsed timesteps * number of atoms * 3 coordinates)
        const LAMMPSFrame* current_frame = data->frames[t];

        for (int p = 0; p < data->num_atoms; p++) {
            const int64_t at_id = data->atomIds[p];
            if(at_id == atom_id) {
                // atom offset = (already considered atoms) * (3 coordinates)
                const int64_t a_offset = 3l*p;
                double x_buff = current_frame->coordinates[a_offset  ];
                double y_buff = current_frame->coordinates[a_offset+1];
                double z_buff = current_frame->coordinates[a_offset+2];
                rarray_push(coordinates_buff, &x_buff);
                rarray_push(coordinates_buff, &y_buff);
                rarray_push(coordinates_buff, &z_buff);
            }
        }
    }

    double* beadTraj = (double*) rarray_to_array(coordinates_buff);
    rarray_free(coordinates_buff);

    return beadTraj;
}



//**************** EXTRACT ATOM IDS FROM MOLECULE AND TYPES ****************
int64_t* getAtomsFromMoleculeList(const LAMMPSData* data,
                                  const int64_t* molecule_ids,
                                  const int64_t numMolecules,
                                  int64_t* numOfSelectedAtoms) {
    /*
     * Return the atom ids corresponding to a given molecule.
     * Supports single molId by passing an array of a single value.
     */
    *numOfSelectedAtoms = 0;
    rarray* atoms_buffer = rarray_init(sizeof(int64_t), 10);

    for(int m=0; m < numMolecules; m++) {
        for(int64_t i=0; i<data->num_atoms; i++) {
            if(data->moleculeIds[i] == molecule_ids[m]) {
                rarray_push(atoms_buffer, &data->atomIds[i]);
                *numOfSelectedAtoms += 1;
            }
        }
    }

    int64_t* sel_atIds = (int64_t*) rarray_to_array(atoms_buffer);
    rarray_free(atoms_buffer);

    return sel_atIds;
}

int64_t* getAtomsFromTypeList(const LAMMPSData* data, const int64_t* a_types, const int64_t numATypes, int64_t* numOfSelectedAtoms) {
    /*
     * Return the atom ids corresponding to a given molecule.
     * Supports single atom type by passing an array of a single value.
     */
    rarray* atoms_buff   = rarray_init(sizeof(int64_t), 10);
    *numOfSelectedAtoms = 0;

    for(int64_t t=0; t<numATypes; t++) {
        for(int64_t i=0; i<data->num_atoms; i++) {
            if(data->atomTypes[i] == a_types[t]) {
                rarray_push(atoms_buff, &data->atomIds[i]);
                *numOfSelectedAtoms += 1;
            }
        }
    }

    if (*numOfSelectedAtoms == 0) {
        rarray_free(atoms_buff);
        return NULL;
    }

    int64_t* sel_atIds = (int64_t*) rarray_to_array(atoms_buff);
    rarray_free(atoms_buff);

    return sel_atIds;
}
