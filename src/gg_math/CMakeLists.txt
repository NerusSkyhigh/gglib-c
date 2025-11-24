# Create a library target from msd.c
add_library(gg_math STATIC
        gg_math.c
)

# Let targets that link to this access its headers
target_include_directories(gg_math
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)