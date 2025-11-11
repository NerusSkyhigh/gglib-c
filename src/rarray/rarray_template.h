// rarray_template.h
/**
 * @brief Template for generating type-specific rarray.
 * Usage: DEFINE_RARRAY(IntArray, int)
 *
 * This macro generates:
 * - struct PingPongArray_NAME
 * - init, free, push, access, size, to_array functions
 *
 * It keeps the same interface as the generic rarray but avoids void* casts.
 */
#define DEFINE_RARRAY(NAME, TYPE) \
typedef struct rarray_##NAME { \
    double growthFactor; \
    size_t numElements; \
    size_t bufferSize; \
    TYPE* data; \
} rarray_##NAME; \
\
static inline rarray_##NAME* rarray_init_##NAME(size_t initialCapacity) { \
    rarray_##NAME* arr = malloc(sizeof(rarray_##NAME)); \
    if (!arr) return NULL; \
    arr->growthFactor = 1; \
    arr->numElements = 0; \
    arr->bufferSize = initialCapacity; \
    arr->data = malloc(sizeof(TYPE) * initialCapacity); \
    if (!arr->data) { free(arr); return NULL; } \
    return arr; \
} \
\
static inline void rarray_free_##NAME(rarray_##NAME* arr) { \
    if (arr) { free(arr->data); free(arr); } \
} \
\
static inline int rarray_push_##NAME(rarray_##NAME* arr, TYPE value) { \
    if (!arr) return -1; \
    if (arr->numElements+1 > arr->bufferSize) { \
        size_t increase = (size_t)((double)arr->bufferSize * arr->growthFactor); \
        size_t newSize = arr->bufferSize + increase; \
        TYPE* newData = realloc(arr->data, sizeof(TYPE) * newSize); \
        if (!newData) return -1; \
        arr->data = newData; \
        arr->bufferSize = newSize; \
    } \
    arr->data[arr->numElements++] = value; \
    return 0; \
} \
\
static inline TYPE* rarray_access_##NAME(rarray_##NAME* arr, size_t index) { \
    if (!arr || index >= arr->numElements) return NULL; \
    return &arr->data[index]; \
} \
\
static inline size_t rarray_size_##NAME(rarray_##NAME* arr) { \
    return arr ? arr->numElements : 0; \
} \
\
static inline TYPE* rarray_to_array_##NAME(rarray_##NAME* arr) { \
    if (!arr) return NULL; \
    TYPE* copy = malloc(sizeof(TYPE) * arr->numElements); \
    if (!copy) return NULL; \
    for (size_t i=0;i<arr->numElements;i++) copy[i] = arr->data[i]; \
    return copy; \
}