# Create a library target from msd.c
add_library(CircularBuffer STATIC
        CircularBuffer.c
)

# Let targets that link to this access its headers
target_include_directories(CircularBuffer
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)