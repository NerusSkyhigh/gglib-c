// CircularBuffer.h
#pragma once
#include <string.h>


/**
 * @struct CircularBuffer
 * @brief Buffer to store the last `capacity` items
 * @param capacity the maximum number of items to store
 * @remarks https://en.wikipedia.org/wiki/Circular_buffer
 */
typedef struct CircularBuffer {
    size_t capacity;
    size_t i;
    size_t count;
    double* data;
} CircularBuffer;

/** @fn initCircularBuffer(CircularBuffer* cb, int64_t capacity)
 * @brief init the memory associateds to the cyclic buffer
 * @param cb pointer to the cyclic buffer to init
 * @param capacity the capacity of the buffer
 */
void initCircularBuffer(CircularBuffer* cb, const size_t capacity);

/** @fn freeCircularBuffer(CircularBuffer* cb)
 * @brief frees the memory associateds to the cyclic buffer
 * @param cb the cyclic buffer to free
 * @warning does NOT free the struct itself
 */
void freeCircularBuffer(CircularBuffer* cb);


/**
 *@brief Push a new element into the buffer. If the buffer is full, the oldest item is replaced
 * @param cb The pointer to the CircularBuffer where to push the element
 * @param e the element to push
 * @warning As the oldest element is replaced, there is not garantee that latest element is the last one
 */
#define CIRCULARBUFFER_PUSH(cb, e) do { \
    (cb)->data[ (cb)->i % (cb)->capacity] = e; \
    (cb)->i++; \
    if( (cb)->count < (cb)->capacity) (cb)->count++; \
} while (0)

/**
 * @brief Returns the last element pushed in the buffer
 * @param cb The pointer to the CircularBuffer where to read the element
 * @warning If no element was stored a `NULL` is returned.
 */
#define CIRCULARBUFFER_LAST(cb) do { \
    if( (cb)->count == 0) return NULL; \
    size_t idx = ( (cb)->i - 1) % (cb)->capacity; \
    return &(cb)->data[idx]; \
} while(0)

/**
 * @brief Returns the oldest element available pushed in the buffer
 * @param cb The pointer to the CircularBuffer where to read the element
 * @warning If no element was stored a `NULL` is returned.
 */
#define  CIRCULARBUFFER_FIRST(cb) do { \
    if ( (cb)->count == 0) return NULL; \
    size_t idx = ( (cb)->i + (cb)->capacity - (cb)->count) % (cb)->capacity; \
    return &(cb)->data[idx]; \
} while(0)




#include "CircularBuffer_template.h"