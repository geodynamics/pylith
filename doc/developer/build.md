# PyLith

## Build developer version with PyLith Installer

### Darwin Mac OS X

```
PATH_TO_INSTALLER_SOURCE_DIR/configure \
    --with-pylith-git=$BRANCH --with-pylith-repo=$REPO_URL \
    --enable-autotools \
	--enable-mpi=mpich \
	--enable-pcre --enable-swig \
	--with-fetch=curl \
    --with-make-threads=2 \
    --with-petsc-options="--download-metis=1 --download-ml=1" \
    --prefix=$HOME/pylith
```

For Poroelasticity group:

    * `BRANCH=josimar/features-poroelasticity`
    * `REPO_URL=https://github.com/josimarsilva/pylith`

## Setup personal fork (do this once per machine)

In the top-level PyLith source directory, run

```
git remote add upstream https://github.com/geodynamics/pylith.git
```

## Update personal fork from geodynamics/pylith

In the top-level PyLith source directory, run

```
MAINBRANCH=knepley/feature-petsc-fe
git fetch upstream
git checkout $MAINBRANCH
git merge upstream/$MAINBRANCH
```

Normally, `MAINBRANCH=master`; however, for current work on the `knepley/feature-petsc-fe` branch.

Merge the main branch into your local branch:

```
git checkout $MYBRANCH
git merge $MAINBRANCH
```

If merge is successful, then push the update to your personal fork:

```
git push
```

## Update local clone of your personal fork on your machine

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

