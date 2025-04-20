cmake_minimum_required(VERSION 3.30)
project(Oscillatore_Smorzato C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -lm -I include") # Match GCC options


add_executable(Oscillatore_Smorzato src/oscill.c)

target_include_directories(Oscillatore_Smorzato PRIVATE include)
target_include_directories(Oscillatore_Smorzato PRIVATE src)