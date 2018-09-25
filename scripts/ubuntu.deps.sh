#!/bin/bash
sudo apt-get install cmake make libgmp-dev python2.7-dev python-sympy python3-dev python3-pip

# No python3-sympy package in Debian
sudo python3 -m pip install sympy
