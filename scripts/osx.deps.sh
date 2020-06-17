#!/bin/bash

echo 'eval "$(pyenv init -)"' >> ${HOME}/.bash_profile

source ${HOME}/.bash_profile

pyenv versions | grep ${PYTHON}
if [[ $? != 0 ]]; then
    pyenv install ${PYTHON}
fi

echo "Using ${PYTHON}"
pyenv global ${PYTHON}

pip install sympy
