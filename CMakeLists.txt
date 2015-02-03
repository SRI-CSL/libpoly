cmake_minimum_required(VERSION 2.6)

project(libpoly C)

# Configure the library source
add_subdirectory(src)

# Configure the headers
add_subdirectory(include)

# Configure the python bindings
add_subdirectory(python)

# Configure the tests
add_subdirectory(test)
