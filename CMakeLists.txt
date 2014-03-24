cmake_minimum_required(VERSION 2.6)

project(libpoly C)

# Default is release with debug info
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF()

# Configure the library source
add_subdirectory(src)

# Conficure the python bindings
add_subdirectory(python)

