// PingPongBuffer_template.h
#pragma once
#include <stdlib.h>
#include <string.h>

/**
 * @brief Template for generating type-specific ping-pong buffers.
 * Usage: DEFINE_PINGPONGBUFFER(NAME, TYPE)
 */
#define DEFINE_PINGPONGBUFFER(NAME, TYPE) \
typedef struct { \
    size_t N; \
    TYPE *data; \
    TYPE *prev; \
    TYPE *next; \
} NAME; \
\
static inline void NAME##_init(NAME *ppb, size_t n_elements) { \
    ppb->N = n_elements; \
    ppb->data = malloc(sizeof(TYPE) * 2 * n_elements); \
    ppb->prev = ppb->data; \
    ppb->next = ppb->data + n_elements; \
} \
\
static inline void NAME##_free(NAME *ppb) { \
    free(ppb->data); \
    ppb->data = NULL; \
    ppb->N = 0; \
} \
\
static inline void NAME##_next_step(NAME *ppb) { \
    TYPE *tmp = ppb->prev; \
    ppb->prev = ppb->next; \
    ppb->next = tmp; \
} \
\
static inline void NAME##_clear_prev(NAME *ppb) { \
    memset(ppb->prev, 0, sizeof(TYPE) * ppb->N); \
} \
\
static inline void NAME##_clear_next(NAME *ppb) { \
    memset(ppb->next, 0, sizeof(TYPE) * ppb->N); \
}
