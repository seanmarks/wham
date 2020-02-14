#!/bin/bash -eE

export CC=gcc
export CXX=g++

# dlib installation
export DLIB_DIR="${HOME}/programs/dlib/19.4"


############################################


### Configure ###

# Create/reset build directory
build_dir="$PWD/build"
if [[ -d $build_dir ]]; then
	rm -r $build_dir
fi
mkdir -p $build_dir

# Configure from build directory
cd $build_dir
cmake .. \
	-DCMAKE_INSTALL_PREFIX="${HOME}/programs/wham" \
	-DCMAKE_PREFIX_PATH="$DLIB_DIR" \

#	-DDLIB_DIR="$DLIB_DIR" \
#	-DOPENMP_ENABLED=ON

### Compile ###

make -j 8


### Install ###

#make install
