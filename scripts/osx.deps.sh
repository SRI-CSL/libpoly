#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp

brew install python2
brew install python3

python2 -m pip install sympy
python3 -m pip install sympy
