#!/bin/bash
brew outdated cmake || brew upgrade cmake
brew outdated gmp || brew upgrade gmp

#
# HACK
#
# Some versions of brew instal python2.7 as python2
# but do not create /usr/local/bin/python.
#
# Because there doesn't seem to be way to force cmake to
# use the python we want (i.e., not /usr/bin/python), we
# create the missing symbolic link here.
#
# This all assumes that /usr/local/bin is before /usr/bin
# on the PATH, which seems to be true for travis.
#
for p in python python-config ;
do
    q=`echo $p | sed -e 's/python/python2/' `
    default=`which $p`
    variant=`which $q`
    if [ "$default" == "/usr/bin/$p" ] && [ "x$variant" != "x" ] ; then
       ln -s $variant /usr/local/bin/$p
    fi
done

#
# Make sure we also have pip
#
default=`which pip`
variant=`which pip2`
if [ "x$default" == "x" ] && [ "x$variant" != "x" ] ; then
    ln -s $variant /usr/local/bin/pip
fi

pip install sympy
