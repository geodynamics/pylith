#!/bin/bash

action="build"
if [ $# == 1 ]; then
  action=$1
fi

if [ $action == "build" ]; then

git clone --depth 1 --recursive https://github.com/geodynamics/pylith_installer.git

pushd pylith_installer && \
  autoreconf --install --verbose && \
  ./configure --prefix=${DEPS_DIR} --enable-mpi=mpich --enable-cppunit --enable-numpy --enable-proj4 --enable-hdf5 --enable-netcdfpy --enable-cmake --enable-setuptools --enable-nemesis --enable-fiat --enable-pcre --enable-swig --with-fortran=no --enable-force-install --with-make-threads=${MAKE_THREADS} && \
  . ./setup.sh && \
  make deps && \
  popd

fi

# End of file
