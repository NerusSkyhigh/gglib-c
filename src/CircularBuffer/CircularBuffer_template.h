#pragma once
#include <stdlib.h>

/**
 * @brief Template macro for generating type-specific circular buffers.
 * Usage: DEFINE_CIRCULARBUFFER(NAME, TYPE)
 *
 * Example:
 *   DEFINE_CIRCULARBUFFER(Double, double)
 *   CircularBuffer_Double cb;
 *   initCircularBuffer_Double(&cb, 32);
 *   CIRCULARBUFFER_PUSH_Double(&cb, 1.23);
 *   printf("%f\n", *CIRCULARBUFFER_LAST_Double(&cb));
 *   freeCircularBuffer_Double(&cb);
 */
#define DEFINE_CIRCULARBUFFER(NAME, TYPE) \
typedef struct CircularBuffer_##NAME { \
    size_t capacity; \
    size_t i; \
    size_t count; \
    TYPE* data; \
} CircularBuffer_##NAME; \
\
static inline void initCircularBuffer_##NAME(CircularBuffer_##NAME* cb, size_t capacity) { \
    cb->capacity = capacity; \
    cb->i = 0; \
    cb->count = 0; \
    cb->data = malloc(sizeof(TYPE) * capacity); \
} \
\
static inline void freeCircularBuffer_##NAME(CircularBuffer_##NAME* cb) { \
    free(cb->data); \
    cb->data = NULL; \
    cb->capacity = 0; \
    cb->count = 0; \
} \
\
static inline void CIRCULARBUFFER_PUSH_##NAME(CircularBuffer_##NAME* cb, TYPE value) { \
    cb->data[cb->i % cb->capacity] = value; \
    cb->i++; \
    if (cb->count < cb->capacity) cb->count++; \
} \
\
static inline TYPE* CIRCULARBUFFER_LAST_##NAME(CircularBuffer_##NAME* cb) { \
    if (cb->count == 0) return NULL; \
    size_t idx = (cb->i - 1) % cb->capacity; \
    return &cb->data[idx]; \
} \
\
static inline TYPE* CIRCULARBUFFER_FIRST_##NAME(CircularBuffer_##NAME* cb) { \
    if (cb->count == 0) return NULL; \
    size_t idx = (cb->i + cb->capacity - cb->count) % cb->capacity; \
    return &cb->data[idx]; \
} \
\
static inline TYPE* CIRCULARBUFFER_GET_##NAME(CircularBuffer_##NAME* cb, size_t offset) { \
    if (offset >= cb->count) return NULL; \
    size_t idx = (cb->i + cb->capacity - cb->count + offset) % cb->capacity; \
    return &cb->data[idx]; \
}
