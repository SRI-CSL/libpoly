#!/bin/bash

pushd .
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar xvf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --host=x86_64-w64-mingw32 --build=--build=$(gcc -dumpmachine) --enable-cxx
make
make install
popd


