// gg_mem.h
// Created by Guglielmo Grillo on 17/10/25.
//
#pragma once
#include <string.h>

/** @struct PingPongBuffer
 * @brief PingPongBuffer
 * @param N Size of a single buffer. Effective length is 2xN
 * @param data All the data in the struct, both old and new values. Total size is 2N
 * @param old Reference to the start of the buffer of old data
 * @param new Reference to the start of the buffer of new data
 * @warning Memory is not cleared on init. Use `PINGPONG_BUFFER_CLEAR_NEXT` and `PINGPONG_BUFFER_CLEAR_PREV`
 * @remark For the pattern see https://gameprogrammingpatterns.com/double-buffer.html
 */
typedef struct PingPongBuffer {
    size_t N;
    double* data;
    double* prev;
    double* next;
} PingPongBuffer;

/**
 *@brief Prepare the pingpong buffer for the next iteration
 * @param db The pointer to the PingPongBuffer to prepare for the next step
 */
#define PINGPONG_BUFFER_NEXT_STEP(ppb) do { \
    double* tmp = (ppb)->prev; \
    (ppb)->prev  = (ppb)->next; \
    (ppb)->next = tmp; \
} while (0)

/**
 *@brief Clear the writing buffer for accumulation in the next step
 * @param db The pointer to the PingPongBuffer where the next buffer is stored
 */
#define PINGPONG_BUFFER_CLEAR_NEXT(ppb) do {\
    memset((ppb)->next, 0, sizeof(double)*(ppb)->N);\
} while (0)

/**
 *@brief Clear the prev buffer
 * @param db The pointer to the PingPongBuffer where the prev buffer is stored
 */
#define PINGPONG_BUFFER_CLEAR_PREV(ppb) do {\
    memset((ppb)->prev, 0, sizeof(double)*(ppb)->N);\
} while (0)

/** @fn initPingPongBuffer(PingPongBuffer* db, int64_t n_elements)
 * @brief init the memory associateds to the ping pong buffer
 * @param ppb the pingpong buffer to init
 * @param n_elements the number of elements for each buffer
 * @warning Memory is not cleared on init. Use `PINGPONG_BUFFER_CLEAR_NEXT` and `PINGPONG_BUFFER_CLEAR_PREV`
 */
void initPingPongBuffer(PingPongBuffer* ppb, const size_t n_elements);

/** @fn freePingPongBuffer(PingPongBuffer* dba)
 * @brief frees the memory associateds to the pingpong buffer
 * @param ppb the pingpong buffer to free
 * @warning does NOT free the struct itself
 */
void freePingPongBuffer(PingPongBuffer* ppb);


#include "PingPongBuffer_template.h"