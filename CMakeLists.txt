cmake_minimum_required(VERSION 2.6)

project(libpoly C)

# Add warnings as errors for good practice
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -pedantic -std=gnu99")

# Configure the library source
add_subdirectory(src)

# Configure the python bindings
add_subdirectory(python)

# Configure the tests
add_subdirectory(test)
