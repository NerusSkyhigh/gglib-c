//
// Created by gu on 14/04/25.
//

#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <stdint.h>
#include <toml.h>

struct config {
    int64_t T_EQ;           // Equilibration time

    int64_t* atom_ids;      // Number of elements in mol_ids
    int64_t  atom_ids_size;

    int64_t* mol_ids;       // List of molecule IDs
    int64_t  mol_ids_size;

    int64_t* atom_types;
    int64_t  atom_types_size;

    char* input_file;
    char* g1_output_file;
    char* g2_output_file;
    char* g3_output_file;
    char* com_output_file;
};

// Function to parse the TOML file and store the values in the config struct
int read_config(const char* filename, struct config* cfg);
void read_int64_array_from_field(toml_table_t *toml_table, const char* field, int64_t** array, int64_t* size);
char* read_filename_from_field(toml_table_t* table, const char* key);
void free_config(struct config* cfg);

#endif // CONFIGPARSER_H

