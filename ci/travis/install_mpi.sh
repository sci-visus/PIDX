#!/bin/bash

# for macOS builds use OpenMPI from homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    cd openmpi
    # check to see if OpenMPI is cached from previous build
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
        echo "Installing OpenMPI with homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/openmpi
        brew install open-mpi
    fi
else

  if [ -f mpich/lib/libmpich.so ]; then
    echo "libmpich.so found -- nothing to build."
  else
    echo "Downloading mpich source."
    wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar xfz mpich-3.2.tar.gz
    rm mpich-3.2.tar.gz
    echo "configuring and building mpich."
    cd mpich-3.2
    ./configure \
        --prefix=`pwd`/../mpich \
        --enable-static=false \
        --enable-alloca=true \
        --disable-long-double \
        --enable-threads=single \
        --enable-fortran=no \
        --enable-fast=all \
        --enable-g=none \
        --enable-timing=none
    make -j4
    make install
    cd -
    rm -rf mpich-3.2
fi

    test -n $CC && unset CC
    test -n $CXX && unset CXX
fi
