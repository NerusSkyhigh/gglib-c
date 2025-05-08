# Create a library target from msd.c
add_library(utils
                memutils.c memutils.h
                fileutils.c fileutils.h
)

# Let targets that link to this access its headers
target_include_directories(utils
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
