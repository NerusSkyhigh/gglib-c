# Create a library target from msd.c
add_library(PingPongBuffer STATIC
        PingPongBuffer.c
)

# Let targets that link to this access its headers
target_include_directories(PingPongBuffer
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)