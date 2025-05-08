#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "rarray.h"


rarray* rarray_init(const size_t item_size, const size_t initialCapacity) {
    if (item_size == 0) {
        fprintf(stderr, "Error: Invalid item size in rarray_init.\n");
        return NULL; // Invalid item size
    } else if (
        initialCapacity == 0) {
        fprintf(stderr, "Error: Invalid initial capacity in rarray_init.\n");
        return NULL; // Invalid initial capacity
    }   

    rarray* arr = malloc(sizeof(rarray));
    if (!arr) {
        fprintf(stderr, "Error: Memory allocation failed in rarray_finalize.\n");
        return NULL; // Allocation failed
    }

    arr->growthFactor = 1; // Default growth factor
    arr->item_size = item_size;
    arr->numElements = 0;
    arr->bufferSize = initialCapacity;
    arr->data = malloc(item_size * initialCapacity);
    
    if (!arr->data) { 
        fprintf(stderr, "Error: Memory allocation failed in rarray_finalize.\n");
        free(arr);
        return NULL;
    }

    return arr;
}


void rarray_free(rarray* arr) {
    if (arr) {
        if (arr->data) {
            free(arr->data);
        }
        free(arr);
    }
}


int rarray_push(rarray* arr, void* newElement) {

    if (!arr || !newElement) {
        fprintf(stderr, "Error: Invalid input in rarray_append.\n");
        return -1;
    }

    // If the element would overflow the buffer, resize the buffer
    if (arr->numElements+1 > arr->bufferSize) {
        if(arr->growthFactor <= 0) {
            fprintf(stderr, "Error: Growth factor must be greater than 0 in rarray_append.\n");
            return -1;
        }

        size_t increase = (size_t) ceil( (double)arr->bufferSize*arr->growthFactor);
        size_t newSize =  (size_t) (arr->bufferSize + increase);
        void*  newData = realloc(arr->data, newSize * arr->item_size);
        if (!newData) {
            fprintf(stderr, "Error: Memory allocation failed in rarray_append.\n");
            return -1;
        }

        arr->data = newData;
        arr->bufferSize = newSize;
    }

    // The cast to char* is necessary to perform pointer arithmetic as (char*) is
    // the only data type guaranteed to be 1 byte in size in C.
    // No need to store the result in a size_t, we can directly use the pointer:
    void* destination = (char*)arr->data + arr->numElements * arr->item_size;
    memcpy(destination, newElement, arr->item_size);

    arr->numElements++;

    return 0;
}

void* rarray_access(rarray* arr, const size_t index) {
    if (!arr) {
        return NULL;
    } else if(index >= arr->numElements) {
        fprintf(stderr, "Error: Index %zu out of bounds in rarray_access.\n", index);
        return NULL;
    }
    return (char*)arr->data + index * arr->item_size; // Return a pointer to the element
}

size_t rarray_size(rarray* arr) {
    if (!arr) {
        fprintf(stderr, "Error: Invalid input in rarray_size.\n");
        return -1;
    }
    return arr->numElements; // Return the number of elements
}



void* rarray_to_array(rarray* arr) {
    if (!arr) {
        fprintf(stderr, "Error: Invalid input in rarray_to_array.\n");
        return NULL;
    } else if (!arr->data) {
        fprintf(stderr, "Error: Invalid data buffer in rarray_to_array.\n");
        return NULL;
    }

    // Allocate new memory to store only the actual elements (not the full buffer)
    void* newArray = malloc(arr->numElements * arr->item_size);
    if (!newArray) {
        fprintf(stderr, "Error: Memory allocation failed in rarray_to_array.\n");
        return NULL;
    }

    memcpy(newArray, arr->data, arr->numElements * arr->item_size);

    return newArray;
}