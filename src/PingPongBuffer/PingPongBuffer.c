// PingPongBuffer.c
// Created by Guglielmo Grillo on 17/10/25.
//
#include <stdlib.h>

#include "PingPongBuffer.h"


void initPingPongBuffer(PingPongBuffer* ppb, const size_t n_elements) {
    ppb->N = n_elements;
    ppb->data = malloc(sizeof(double)*2*n_elements);
    ppb->prev = &(ppb->data[0]);
    ppb->next = &(ppb->data[n_elements]);
}

void freePingPongBuffer(PingPongBuffer* ppb) {
    free(ppb->data);
    ppb->N = 0;
}
