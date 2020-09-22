#!/bin/bash

pushd .
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar xvf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --host=x86_64-w64-mingw32 --enable-cxx --prefix=$HOME/prefix
make
make install
popd


