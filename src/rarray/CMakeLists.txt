# Create a library target from msd.c
add_library(rarray rarray.c)

# Let targets that link to this access its headers
target_include_directories(rarray
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
