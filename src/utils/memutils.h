// memutils.h
#pragma once
#include <stdio.h>
#include <stdlib.h>

#define SAFE_MALLOC(ptr, size, msg) do {                                                        \
    ptr = malloc(size);                                                                         \
    if (!(ptr)) {                                                                               \
        fprintf(stderr, "[%s:%d] %s: malloc failed; %s\n", __FILE__, __LINE__, __func__, msg);  \
        exit(EXIT_FAILURE);                                                                     \
    }                                                                                           \
} while (0)

#define SAFE_CALLOC(ptr, count, size, msg) do {                                                 \
    ptr = calloc(count, size);                                                                  \
    if (!(ptr)) {                                                                               \
        fprintf(stderr, "[%s:%d] %s: calloc failed; %s\n", __FILE__, __LINE__, __func__, msg);  \
        exit(EXIT_FAILURE);                                                                     \
    }                                                                                           \
} while (0)

#define SAFE_FREE(ptr) do {     \
    if (ptr) {                  \
        free(ptr);              \
        ptr = NULL;             \
    }                           \
} while (0)

#define SAFE_REALLOC(ptr, size, msg) do {                                                       \
    void* tmp = realloc(ptr, size);                                                             \
    if (!tmp) {                                                                                 \
        fprintf(stderr, "[%s:%d] %s: realloc failed; %s\n", __FILE__, __LINE__, __func__, msg); \
        exit(EXIT_FAILURE);                                                                     \
    }                                                                                           \
    ptr = tmp;                                                                                  \
} while (0)
