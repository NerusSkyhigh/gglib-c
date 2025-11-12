// CircularBuffer.c
// Created by Gugli on 12/11/2025.
//

#include "CircularBuffer.h"
#include <stdlib.h>


void initCircularBuffer(CircularBuffer* cb, const size_t capacity) {
    cb->data = malloc(sizeof(double)*capacity);
    cb->capacity = capacity;    // [TODO] Safeguard allocation
    cb->i = 0;
    cb->count = 0;
}

void freeCircularBuffer(CircularBuffer* cb) {
    free(cb->data);
    cb->capacity = 0;
    cb->count = 0;
}