// rarray.h
#pragma once
#include <stddef.h>

/**
 * @struct rarray
 * @brief Generic runtime-resizable array using void pointers.
 * @warning The default rarray uses void*, so you must cast elements when reading/writing.
 */
typedef struct rarray {
    double growthFactor; ///< Growth factor for resizing the buffer
    size_t numElements;  ///< Number of elements in the array
    size_t item_size;    ///< Size of each element in bytes
    size_t bufferSize;   ///< Allocated buffer size in elements
    void* data;          ///< Pointer to the data buffer
} rarray;

/**
 * @brief Initialize a generic rarray
 * @param item_size Size of each element in bytes
 * @param initialCapacity Initial number of elements to allocate
 * @return Pointer to the allocated rarray
 */
rarray* rarray_init(size_t item_size, size_t initialCapacity);

/**
 * @brief Free memory allocated for rarray
 * @param arr Pointer to the rarray
 */
void rarray_free(rarray* arr);

/**
 * @brief Append a new element to the rarray
 * @param arr Pointer to the rarray
 * @param newElement Pointer to the element to append
 * @return 0 on success, -1 on failure
 */
int rarray_push(rarray* arr, void* newElement);

/**
 * @brief Access an element of the rarray
 * @param arr Pointer to the rarray
 * @param index Index of the element
 * @return Pointer to the element, NULL if out of bounds
 */
void* rarray_access(rarray* arr, size_t index);

/**
 * @brief Return the number of elements in the rarray
 * @param arr Pointer to the rarray
 * @return Number of elements
 */
size_t rarray_size(rarray* arr);

/**
 * @brief Copy the active elements to a new array
 * @param arr Pointer to the rarray
 * @return Pointer to newly allocated array containing all elements
 */
void* rarray_to_array(rarray* arr);

#include "rarray_template.h"