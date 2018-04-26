#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp
brew install python
pip2 install sympy