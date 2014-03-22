cmake_minimum_required(VERSION 2.6)

project(libpoly C)

# Default is release with debug info
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF()

# Find the GMP number library
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
if (GMP_INCLUDE_DIR AND GMP_LIBRARY)
    message(STATUS "GMP headers: "${GMP_INCLUDE_DIR})
    message(STATUS "GMP library: "${GMP_LIBRARY})
    include_directories(${GMP_INCLUDE_DIR})
else()
    MESSAGE(FATAL_ERROR "Could not the GMP number library (sudo apt-get install libgmp-dev)")
endif()

# Add our sources to the include
include_directories(${libpoly_SOURCE_DIR}/src ${libpoly_SOURCE_DIR}/include)

# Configure the subdirectories
add_subdirectory(src)

