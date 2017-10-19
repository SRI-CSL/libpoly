#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp
brew outdated python || brew upgrade python
sudo easy_install pip
pip install sympy