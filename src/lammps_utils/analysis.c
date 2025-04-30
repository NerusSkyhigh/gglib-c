//
// Created by gu on 27/03/25.
//
#include <stdint.h>
#include <stdio.h>

#include "analysis.h"
#include "parser.h"
#include "rarray.h"

void compute_CoM(const double* frame, const int64_t num_atoms, double* com_coord) {
    // Set to zero for accumulation
    com_coord[0] = 0.;
    com_coord[1] = 0.;
    com_coord[2] = 0.;

    for (int i = 0; i < num_atoms; i++) { // TODO: Transform this into a `size_t` object
        // Let's divide by number of atoms to avoid overflow
        com_coord[0] += frame[3*i]   / ( (double) num_atoms);
        com_coord[1] += frame[3*i+1] / ( (double) num_atoms);
        com_coord[2] += frame[3*i+2] / ( (double) num_atoms);
    }
    //printf("%lf %lf %lf\n", com_coord[0], com_coord[1], com_coord[2]);
}

double* getBeadTrajectory(const LAMMPSData* data, const int64_t atom_id) {
    /*
     * Return the trajectory of a select bead as a function of time.
     * The array is a linear array.
     */
    rarray* coordinates_buff   = rarray_init(sizeof(double), 10);

    for (int64_t t = 0; t < data->num_timesteps; t++)
    {
        // Frame offset (elapsed timesteps * number of atoms * 3 coordinates)
        const int64_t f_offset = t * data->num_atoms * 3;

        for (int p = 0; p < data->num_atoms; p++) {
            const int64_t at_id = data->atomIds[p];

            if(at_id == atom_id) {
                // atom offset = (already considered atoms) * (3 coordinates)
                const int64_t a_offset = 3*p;
                double x_buff = data->coordinates[f_offset+a_offset  ];
                double y_buff = data->coordinates[f_offset+a_offset+1];
                double z_buff = data->coordinates[f_offset+a_offset+2];
                rarray_push(coordinates_buff, &x_buff);
                rarray_push(coordinates_buff, &y_buff);
                rarray_push(coordinates_buff, &z_buff);
            }
        }
    }
    //[TODO] Understand how struct are passed to functions. I might be using struct* and struct->value and struct.value badly
    double* beadTraj = (double*) rarray_to_array(coordinates_buff);
    rarray_free(coordinates_buff);

    return beadTraj;
}


int64_t* getAtomsFromMoleculeList(const LAMMPSData* data, const int64_t* molecule_ids, const int64_t numMolecules, int64_t* numOfSelectedAtoms) {

    *numOfSelectedAtoms = 0;
    // We espect at least 1 atom per molecule
    rarray* atoms_buffer = rarray_init(sizeof(int64_t), numMolecules);

    for(int m=0; m < numMolecules; m++) {
        const int64_t molId = molecule_ids[m];
        int64_t numAtomsInMolecule;
        //[TODO] Fix all the pointers to struct around. I'm using `&(*data)` instead of simply `data` as a reminder.
        int64_t* atoms_in_molecule = getAtomsInMolecule(&(*data), molId, &numAtomsInMolecule);

        for(int p=0; p<numAtomsInMolecule; p++) {
            // atoms_in_molecule+p = &(atoms_in_molecule[p])
            rarray_push(atoms_buffer, atoms_in_molecule+p);
            *numOfSelectedAtoms += 1;
        }
        free(atoms_in_molecule);
    }

    int64_t* sel_atIds = (int64_t*) rarray_to_array(atoms_buffer);
    rarray_free(atoms_buffer);

    return sel_atIds;
}

int64_t* getAtomsInMolecule(const LAMMPSData* data, const int64_t molecule_id, int64_t* numOfSelectedAtoms) {
    /*
     * Return the atom ids corresponding to a given molecule.
     */
    rarray* atoms_buff   = rarray_init(sizeof(int64_t), 10);
    *numOfSelectedAtoms = 0;

    for(int64_t i=0; i<data->num_atoms; i++) {
        if(data->moleculeIds[i] == molecule_id) {
            rarray_push(atoms_buff, &data->atomIds[i]);
            *numOfSelectedAtoms += 1;
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

int64_t* getAtomsFromType(const LAMMPSData* data, const int64_t* a_types, const int64_t numATypes, int64_t* numOfSelectedAtoms) {
    /*
     * Return the atom ids corresponding to a given molecule.
     */
    rarray* atoms_buff   = rarray_init(sizeof(int64_t), 10);
    *numOfSelectedAtoms = 0;

    for(int64_t t=0; t<numATypes; t++) {
        int64_t a_type = a_types[t];
        for(int64_t i=0; i<data->num_atoms; i++) {
            if(data->atomTypes[i] == a_type) {
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
