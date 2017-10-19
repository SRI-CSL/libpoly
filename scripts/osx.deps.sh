#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp
virtualenv env
source env/bin/activate
pip install sympy
