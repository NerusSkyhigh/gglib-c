//
// Created by gu on 27/03/25.
//
#ifndef LAMMPSDATAANALYSIS_H
#define LAMMPSDATAANALYSIS_H

#include <stdint.h>

#include "LAMMPSParser.h"

void compute_CoM(const double* frame, const int64_t num_atoms, double* com_coord);
double* getBeadTrajectory(const LAMMPSData* data, const int64_t atom_id);
int64_t* getAtomsFromMoleculeList(const LAMMPSData* data, const int64_t* molecule_ids, const int64_t numMolecules, int64_t* numOfSelectedAtoms);
int64_t* getAtomsFromTypeList(const LAMMPSData* data, const int64_t* a_types, const int64_t numATypes, int64_t* numOfSelectedAtoms);
#endif //LAMMPSDATAANALYSIS_H
