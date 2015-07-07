cmake_minimum_required(VERSION 2.6)

project(libpoly C)

if(APPLE)
 # This is the simplest setting that seems to currently work on the mac.
 # 
 set(CMAKE_MACOSX_RPATH OFF)
 set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
endif()

include(CheckFunctionExists)

check_function_exists(open_memstream HAVE_OPEN_MEMSTREAM)

# Configure the library source
add_subdirectory(src)

# Configure the headers
add_subdirectory(include)

if(HAVE_OPEN_MEMSTREAM)

  # Configure the python bindings	
  add_subdirectory(python)

  # Configure the tests
  add_subdirectory(test)

endif()
