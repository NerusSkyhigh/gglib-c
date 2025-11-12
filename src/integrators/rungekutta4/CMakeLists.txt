cmake_minimum_required(VERSION 3.30)
project(simulator C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

#LIBRARY
add_library(simulator STATIC
        integrators.h
        integrators.c
)
target_include_directories(simulator
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}
)


# ------- Raylib cmake -------
#src: https://github.com/SasLuca/raylib-cmake-template
project(my_raylib_game C)
set(CMAKE_C_STANDARD 99)

#  Adding Raylib
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples
set(BUILD_GAMES    OFF CACHE BOOL "" FORCE) # don't build the supplied example games

FetchContent_Declare(
        raylib
        GIT_REPOSITORY "https://github.com/raysan5/raylib.git"
        GIT_TAG "master"
        GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(raylib)

# Adding our source files
#file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/sources/*.c") # Define PROJECT_SOURCES as a list of all source files
#set(PROJECT_INCLUDE "${CMAKE_CURRENT_LIST_DIR}/sources/") # Define PROJECT_INCLUDE to be the path to the include directory of the project

# Declaring our executable
#add_executable(${PROJECT_NAME})
#target_sources(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCES})
#target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_INCLUDE})
#target_link_libraries(${PROJECT_NAME} PRIVATE raylib)
# ------- END RAYLIB CMAKE -------




# RK4 test case
add_executable(test_RK4
        tests/RK4.c
)

target_link_libraries(test_RK4
        simulator
        raylib
)
IF (NOT WIN32)
    target_link_libraries(test_RK4 m)
ENDIF()