cmake_minimum_required(VERSION 2.6)


project(libpoly C)

include(CheckFunctionExists)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

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
