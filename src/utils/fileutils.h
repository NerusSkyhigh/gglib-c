// fileutils.h
#pragma once
#include <stdio.h>
#include <stdlib.h>

#define SAFE_FILE_OPEN(ptr, filename, modes, msg) do {                      \
    ptr = fopen(filename, modes);                                           \
    if (!(ptr)) {                                                           \
        fprintf(stderr, "[%s:%d] %s: fopen(\"%s\", \"%s\") failed: %s\n",   \
                __FILE__, __LINE__, __func__, filename, modes, msg);        \
        exit(EXIT_FAILURE);                                                 \
    }                                                                       \
} while (0)

#define SAFE_FILE_CLOSE(ptr) do {                               \
    if ((ptr) && fclose(ptr) != 0) {                            \
        fprintf(stderr, "[%s:%d] %s: fclose() failed on %s\n",  \
                __FILE__, __LINE__, __func__, #ptr);            \
        exit(EXIT_FAILURE);                                     \
    }                                                           \
    ptr = NULL;                                                 \
} while (0)

