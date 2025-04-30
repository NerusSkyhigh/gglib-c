//
// Created by gu on 14/04/25.
//

#include "configParser.h"
#include <toml.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int read_config(const char* filename, struct config* cfg) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        perror("Failed to open TOML file");
        return -1;
    }

    // Create a TOML parser
    char errbuf[200];
    toml_table_t *toml_table = toml_parse_file(f, errbuf, sizeof(errbuf));
    fclose(f);

    if (!toml_table) {
        fprintf(stderr, "Failed to parse TOML file: %s\n", filename);
        return -1;
    }

    // READ `T_EQ` (int64_t, required)
    toml_datum_t t_eq = toml_int_in(toml_table, "T_EQ");
    if (!t_eq.ok) {
        fprintf(stderr, "Missing required field: T_EQ\n");
        toml_free(toml_table);
        return -3;
    }
    cfg->T_EQ = t_eq.u.i;

    read_int64_array_from_field(toml_table, "atom_ids",   &cfg->atom_ids,   &cfg->atom_ids_size);
    read_int64_array_from_field(toml_table, "mol_ids",    &cfg->mol_ids,    &cfg->mol_ids_size);
    read_int64_array_from_field(toml_table, "atom_types", &cfg->atom_types, &cfg->atom_types_size);

    // Read filenames from config file.
    // If the field doesn't exist, it will return NULL
    cfg->input_file       = read_filename_from_field(toml_table, "input_file");
    cfg->g1_output_file   = read_filename_from_field(toml_table, "g1_output_file");
    cfg->g2_output_file   = read_filename_from_field(toml_table, "g2_output_file");
    cfg->g3_output_file   = read_filename_from_field(toml_table, "g3_output_file");
    cfg->com_output_file  = read_filename_from_field(toml_table, "com_output_file");

    toml_free(toml_table);
    return 0;
}

void read_int64_array_from_field(toml_table_t *toml_table, const char* field, int64_t** array, int64_t* size) {
    // Init value to 0 as a defaul
    *size = 0;
    *array = NULL;

    // Check if the key exists in the TOML table
    toml_array_t* arr = toml_array_in(toml_table, field);
    if(!arr) return;

    *size = (int64_t) toml_array_nelem(arr);
    *array = malloc(*size * sizeof(int64_t));
    if (!*array) {
        *size = 0;
        fprintf(stderr, "Failed to allocate memory for %s\n", field);
        return;
    }

    for (int64_t i = 0; i < *size; i++) {
        toml_datum_t toml_data = toml_int_at(arr, i);

        if (toml_data.ok) {
             (*array)[i] = (int64_t) toml_data.u.i;
        } else {
            fprintf(stderr, "Warning: could not parse '%s[%ld]'\n", field, i);
            exit(-1);
        }
    }
}

char* read_filename_from_field(toml_table_t* table, const char* key) {
    toml_datum_t val = toml_string_in(table, key);
    if (!val.ok) return NULL;

    size_t len = strlen(val.u.s);
    char* copy = malloc(len + 1);
    if (!copy) {
        fprintf(stderr, "Failed to allocate memory for '%s'\n", key);
        free(val.u.s);
        return NULL;
    }

    memcpy(copy, val.u.s, len + 1);  // includes null terminator
    free(val.u.s);  // must free what toml_string_in allocates
    return copy;
}

void free_config(struct config* cfg) {
    free(cfg->atom_ids);
    free(cfg->mol_ids);
    free(cfg->atom_types);

    free(cfg->input_file);
    free(cfg->g1_output_file);
    free(cfg->g2_output_file);
    free(cfg->g3_output_file);
    free(cfg->com_output_file);
}
