# Create a library target from the lammps_utils folder source files
add_library(lammps_utils
        analysis.c
        parser.c
)

# Specify the directory for the header files
target_include_directories(lammps_utils
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}  # This makes the headers accessible for other targets
)

# Link necessary libraries (if any) here
target_link_libraries(lammps_utils rarray)

# Compiler options can be inherited from the top-level CMake configuration
