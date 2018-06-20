# PyLith

## Update local clone of on your machine

```
git pull
```

If `configure.ac` changes, then run `autoreconf -if` in the top-level PyLith source directory.

## Configure

### With debugging

From the PyLith build directory, run configure:

```
PREFIX=PATH_TO_PYLITH_DEST
PATH_TO_PYLITH_SOURCE/configure --enable-swig --enable-testing --enable-hdf5 --enable-cubit --prefix=$PREFIX CPPFLAGS=-I$PREFIX LDFLAGS=-L$PREFIX CC=mpicc CXX=mpicxx CFLAGS="-g -Wall" CXXFLAGS="-g -Wall"
```

## Build

```
make
make install
```
