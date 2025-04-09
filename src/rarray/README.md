# rarray - A Simple Resizable Array Library

This is a simple C library that implements a resizable array. It provides basic functions to initialize, append, access, and finalize the array.

## Features
- Dynamic resizing with a configurable growth factor.
- Type-agnostic array using `void*`.
- Memory-efficient with proper handling of dynamic memory.

## Installation

### Option 1: Include the Source Files

1. Copy the `include/rarray.h` file to your project's `include/` directory.
2. Copy the `src/rarray.c` file to your project's `src/` directory.
3. Include `rarray.h` in your code:
   ```c
   #include "rarray.h"
