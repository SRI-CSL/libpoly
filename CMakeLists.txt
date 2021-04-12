cmake_minimum_required(VERSION 2.6)

#
# Stop warning on cygwin
#
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

project(libpoly C CXX)

option(LIBPOLY_BUILD_PYTHON_API "Build the Python API" ON)
option(LIBPOLY_BUILD_STATIC_PIC "Build the static PIC library" ON)
option(LIBPOLY_BUILD_STATIC "Build the static library" ON)
option(LIBPOLY_BUILD_STATISTICS "Build the statistics internals" OFF)

set(LIBPOLY_VERSION_MAJOR 0)
set(LIBPOLY_VERSION_MINOR 1)
set(LIBPOLY_VERSION_PATCH 9)
set(LIBPOLY_VERSION ${LIBPOLY_VERSION_MAJOR}.${LIBPOLY_VERSION_MINOR}.${LIBPOLY_VERSION_PATCH})

if(APPLE)
 # This is the simplest setting that seems to currently work on the mac.
 #
 set(CMAKE_MACOSX_RPATH OFF)
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#
# Search for GMP. We do this here since GMP_INCLUDE_DIR may be
# required for building the python binding.
#
# Find the GMP number library
find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
if (GMP_INCLUDE_DIR AND GMP_LIBRARY)
    message(STATUS "GMP headers: " ${GMP_INCLUDE_DIR})
    message(STATUS "GMP library: " ${GMP_LIBRARY})
else()
    MESSAGE(FATAL_ERROR "Could not the GMP number library (sudo apt-get install libgmp-dev)")
endif()


# statistics configuration
if(LIBPOLY_BUILD_STATISTICS)

  add_definitions(-DLIBPOLY_STATISTICS)

endif()

#
# Check for open_memstream
#
include(CheckFunctionExists)

check_function_exists(open_memstream HAVE_OPEN_MEMSTREAM)

if(HAVE_OPEN_MEMSTREAM)
  add_definitions(-DHAVE_OPEN_MEMSTREAM)
endif()

# Configure the library source
add_subdirectory(src)

# Configure the headers
add_subdirectory(include)

if(LIBPOLY_BUILD_PYTHON_API)

  # Need the Python binary to run tests
  find_package(PythonInterp)

  # Configure the Python bindings
  add_subdirectory(python)

  # Configure the Python tests
  add_subdirectory(test/python)

endif()

# Configure the C++ tests
enable_testing()
add_subdirectory(test/polyxx)

