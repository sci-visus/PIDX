#!/bin/bash
##
## BSD 3-Clause License
## 
## Copyright (c) 2010-2018 ViSUS L.L.C., 
## Scientific Computing and Imaging Institute of the University of Utah
## 
## ViSUS L.L.C., 50 W. Broadway, Ste. 300, 84101-2044 Salt Lake City, UT
## University of Utah, 72 S Central Campus Dr, Room 3750, 84112 Salt Lake City, UT
##  
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## * Redistributions of source code must retain the above copyright notice, this
## list of conditions and the following disclaimer.
## 
## * Redistributions in binary form must reproduce the above copyright notice,
## this list of conditions and the following disclaimer in the documentation
## and/or other materials provided with the distribution.
## 
## * Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## For additional information about this project contact: pascucci@acm.org
## For support: support@visus.net
## 
##

# for macOS builds use OpenMPI from homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    cd openmpi
    # check to see if OpenMPI is cached from previous build
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
        echo "Installing OpenMPI with homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/openmpi
        brew install -v open-mpi 
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
