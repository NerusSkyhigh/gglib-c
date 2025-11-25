# Create a library target from the lammps_utils folder source files
add_library(lammpstrjIO STATIC
        lammpstrjIO.c
)

# Specify the directory for the header files
target_include_directories(lammpstrjIO
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
