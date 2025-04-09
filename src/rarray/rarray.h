#ifndef RARRAY_H
#define RARRAY_H

#include <stddef.h>  // For size_t

// The rarray struct definition
struct rarray {
    double growthFactor; // Growth factor for resizing the buffer
    size_t numElements;  // Number of elements in the array
    size_t item_size;    // Size of each element in bytes
    size_t bufferSize;   // Size of the buffer in elements
    void* data;          // Pointer to the data buffer
};
typedef struct rarray rarray;


rarray* rarray_init(size_t item_size, size_t initialCapacity);
void rarray_free(rarray* arr);
int rarray_push(rarray* arr, void* newElement);
void* rarray_access(rarray* arr, size_t index);
size_t rarray_size(rarray* arr);
void* rarray_to_array(rarray* arr);

#endif // RARRAY_H
