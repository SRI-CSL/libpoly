#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp
pip2 install sympy