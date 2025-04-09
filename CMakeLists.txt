cmake_minimum_required(VERSION 3.20)
project(BANGLES C)

set(CMAKE_C_STANDARD 11)

# Enable CTest for automated tests
enable_testing()

# Find OpenMP
find_package(OpenMP REQUIRED)


# ──────── Subdirectories ────────
add_subdirectory(src)
add_subdirectory(tests)


# ──────── Main executable ────────
add_executable(BANGLES computeMSD.c)


# ──────── Link with modular libraries ────────
target_link_libraries(BANGLES
        PRIVATE
        msd
        lammps_utils
        rarray
        OpenMP::OpenMP_C
        m                 # math library
)


# ──────── Keep warnings strict ────────
if(MSVC)
    target_compile_options(BANGLES PRIVATE /W4 /WX)
else()
    target_compile_options(BANGLES PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()