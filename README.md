# LibPoly

LibPoly is a C library for manipulating polynomials. The target applications 
are symbolic reasoning engines, such as SMT solvers, that need to reason about
polynomial constraints. It is research software under development, so the 
features and the API might change rapidly.

## Prerequisites

To compile on an Ubuntu machine you can install the prerequisites with  

```
sudo apt-get install gcc cmake make libgmp-dev python2.7-dev 
```

Python is used for testing the library purposes through the Python bindings 
(these are also useful for playing with the library). The Python bindings are 
currently not supported for Macs or Windows machines.

## Compiling

To compile and install perform 
```
cd build
cmake .. -DCMAKE_BUILD_TYPE=$type -DCMAKE_INSTALL_PREFIX=$prefix
make
make install
```
The $type above is should be either "Debug" or "Release", and $prefix is the 
target directory where the library will be installed. The prefix can be 
omitted, in which case the library will be installed in the default system 
location (such as ```/usr/local```).

If the tests are enabled, you can do a sanity check of the library by doing a
```make check```.