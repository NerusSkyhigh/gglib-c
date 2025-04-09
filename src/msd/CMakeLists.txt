# Create a library target from msd.c
add_library(msd msd.c)

# Let targets that link to this access its headers
target_include_directories(msd
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(msd rarray)