# src/config_parser/CMakeLists.txt

set(TOML_DIR ${CMAKE_CURRENT_BINARY_DIR}/tomlc99)
file(MAKE_DIRECTORY ${TOML_DIR})

# Download the header and source if they don't already exist
# This is an hack because CMake doesn't work well with external non-cmake libraries.
# I have cmake and my life choices.
file(DOWNLOAD
        https://raw.githubusercontent.com/cktan/tomlc99/master/toml.h
        ${TOML_DIR}/toml.h
        SHOW_PROGRESS
)

file(DOWNLOAD
        https://raw.githubusercontent.com/cktan/tomlc99/master/toml.c
        ${TOML_DIR}/toml.c
        SHOW_PROGRESS
)

# Add tomlc99 directly to your project
add_library(tomlc99 STATIC ${TOML_DIR}/toml.c)
target_include_directories(tomlc99 PUBLIC ${TOML_DIR})

# Your configParser library
add_library(configParser configParser.c)
target_include_directories(configParser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(configParser PRIVATE tomlc99)
