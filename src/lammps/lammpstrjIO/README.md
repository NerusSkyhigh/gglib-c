# LAMMPS TRaJectory Input Output — Minimal LAMMPS Trajectory Writer

This lightweight C library provides functions to create and append frames to a LAMMPS-compatible trajectory file, suitable for visualization in tools such as VMD.
It abstracts away file formatting details, letting you focus on simulation logic.


## Features

- Write trajectories directly in LAMMPS “dump” format (.lammpstrj).
- Automatically handles timestep numbering and box bounds.
- Simple API: initialize → write frames → free.

## 
**Installation**

Just add the two files to your project:
```
lammpstrjIO.c
lammpstrjIO.h
```
or add the project folder to your `CMakeLists.txt`cmake via `add_subdirectory(lammpstrjIO)`.

Then, include the header: `#include "lammpstrjIO.h"` and compile it together with your code: `gcc main.c lammpstrjIO.c -o my_sim`

## API Overview
Check the header files for more details.
- `void initLammpsTrjData(LammpsTrjFile *ld, const char *filename, size_t n_particles, double boxL)`  
Initializes a trajectory descriptor.

- `void writeLammpsTrjFrame(LammpsTrjFile *ld, float *coordinates)`  
Appends a new frame to the file. Automatically creates the file on first call. Coordinates are expected in the format `{x0, y0, z0, x1, y1, z1, ..., xN-1, yN-1, zN-1}` (linear array particle major).

- `void freeLammpsTrjData(LammpsTrjFile *ld)`  
Frees internal memory (e.g., the filename string) and resets counters.

### Example Usage
```c
#include <stdio.h>
#include "gg_io.h"

int main() {
    const size_t N = 4;
    double boxL = 10.0;
    float coords[3*N] = {
        0,0,0,  1,0,0,  0,1,0,  0,0,1
    };
    
    LammpsTrjFile traj;
    initLammpsTrjData(&traj, "out.lammpstrj", N, boxL);

    for (int step = 0; step < 3; step++) {
        for (size_t i = 0; i < N; i++) coords[3*i+0] += 0.1f; // simple motion
        writeLammpsTrjFrame(&traj, coords);
    }

    freeLammpsTrjData(&traj);
    return 0;
}
```

This will create a `out.lammpstrj` file you can open in VMD:

```
Output Format Example
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
4
ITEM: BOX BOUNDS pp pp pp
-5.000000e+00 5.000000e+00
-5.000000e+00 5.000000e+00
-5.000000e+00 5.000000e+00
ITEM: ATOMS id mol type xu yu zu
1 1 1 0.000000 0.000000 0.000000
2 1 1 1.000000 0.000000 0.000000
3 1 1 0.000000 1.000000 0.000000
4 1 1 0.000000 0.000000 1.000000
```

## Notes
- The simulation box is cubic, extending from -boxL/2 to +boxL/2 in each direction.
- All atoms are assigned mol=1, type=1 by default (edit source if you need more complex typing).
- Compatible with standard LAMMPS and VMD readers.