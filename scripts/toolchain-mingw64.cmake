#
# For cross-compilation using cmake.
# This uses the mingw64 tools to build 64bit Windows binaries
# from cygwin. To use this, call cmake with option
#  -DCMAKE_TOOLCHAIN_FILE=<this-file>
#
# 
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_AR /usr/bin/x86_64-w64-mingw32-ar.exe)
set(CMAKE_RANLIB /usr/bin/x86_64-w64-mingw32-ranlib.exe)
set(CMAKE_C_COMPILER /usr/bin/x86_64-w64-mingw32-gcc.exe)
set(CMAKE_CXX_COMPILER /usr/bin/x86_64-w64-mingw32-g++.exe)
set(CMAKE_FIND_ROOT_PATH /home/bruno/tools/mingw64)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
add_definitions(-D__USE_MINGW_ANSI_STDIO)
