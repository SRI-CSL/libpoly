#!/bin/bash
sudo apt-get install cmake make libgmp-dev

echo $PYTHON

echo 'eval "$(pyenv init -)"' >> ${HOME}/.bash_profile;
source ${HOME}/.bash_profile;

# iam: this should show us what our choices are.
pyenv versions

pyenv global ${PYTHON};

pip install sympy
