[![Build Status](https://travis-ci.org/SRI-CSL/libpoly.svg?branch=master)](https://travis-ci.org/SRI-CSL/libpoly)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/5716/badge.svg)](https://scan.coverity.com/projects/5716)


# LibPoly

LibPoly is a C library for manipulating polynomials. The target applications 
are symbolic reasoning engines, such as SMT solvers, that need to reason about
polynomial constraints. It is research software under development, so the 
features and the API might change rapidly.

## Prerequisites

To compile on an Ubuntu machine you can install the prerequisites with  

```
sudo apt-get install gcc cmake make libgmp-dev python2.7-dev python-sympy
```

Python is used for testing the library through the Python bindings 
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


## Installing Prebuilt Binaries

Currently you can install the libpoly library (without python support) either using Homebrew or Aptitude.

#### Homebrew

Installing on Darwin using homebrew can be achieved via:
```
brew install SRI-CSL/sri-csl/libpoly
```


#### Aptitude

To install libpoly on Ubuntu or Debian, do the following:
```
sudo add-apt-repository ppa:sri-csl/formal-methods
sudo apt-get update
sudo apt-get install libpoly-dev
```
This will install the full mcsat enabled version of Yices (as well as the headers and dynamic library).



