#!/bin/bash

action="build"
if [ $# == 1 ]; then
  action=$1
fi

if [ $action == "build" ]; then

git clone --depth 1 --recursive https://github.com/geodynamics/pylith_installer.git

pushd pylith_installer && \
  autoreconf --install --verbose && \
  ./configure --prefix=${DEPS_DIR} --disable-cmake --disable-cppunit --disable-proj4 --enable-numpy --enable-hdf5 --enable-netcdf --enable-netcdfpy --enable-nemesis --enable-fiat --with-fortran=no --enable-force-install --with-make-threads=${MAKE_THREADS} CFLAGS="-O2 -w" CXXFLAGS="-O2 -w" && \
  . ./setup.sh && \
  make external_deps && \
  popd

fi

# End of file
