# Create a library target from the lammps_utils folder source files
add_library(lammpstryIO
        lammpstrjIO.c
)

# Specify the directory for the header files
target_include_directories(lammpstryIO
        PUBLIC
        ..  # This makes the headers accessible for other targets
)