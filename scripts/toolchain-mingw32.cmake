#
# For cross-compilation using cmake.
# This uses the mingw32 tools to build 32bit Windows binaries
# from cygwin. To use this, call cmake with option
#  -DCMAKE_TOOLCHAIN_FILE=<this-file>
#
# Note: on cygwin, there are two mingw tool chains, with
# tools like:
#  /usr/bin/i686-w64-mingw32-gcc.exe
#  /usr/bin/i686-pc-ming32-gcc.exe
#
# These are different and libpoly requires the w64 tool chain.
# 
set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_AR /usr/bin/i686-w64-mingw32-ar.exe)
set(CMAKE_RANLIB /usr/bin/i686-w64-mingw32-ranlib.exe)
set(CMAKE_C_COMPILER /usr/bin/i686-w64-mingw32-gcc.exe)
set(CMAKE_CXX_COMPILER /usr/bin/i686-w64-mingw32-g++.exe)
set(CMAKE_FIND_ROOT_PATH /home/bruno/tools/mingw32)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
add_definitions(-D__USE_MINGW_ANSI_STDIO)
