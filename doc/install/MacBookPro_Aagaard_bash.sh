# Environment variables for Bash Shells associated with Brad Aagaard's
# install of PyLith on a MacBookPro running OS X 10.4

# General
export TOOLS_DIR=/sw/tools
export TOOLS_FORMAT=gcc-4.0
export PYTHON_VERSION=2.5

# Location for installation of CIG software
export CIG_DIR=$HOME/tools/cig/${TOOLS_FORMAT}
export CIG_INCDIR=${CIG_DIR}/include
export CIG_LIBDIR=${CIG_DIR}/lib
PATH=$PATH:${CIG_DIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CIG_LIBDIR}
export PYTHONPATH=${PYTHONPATH}:${CIG_LIBDIR}/python${PYTHON_VERSION}/site-packages

# Python 2.5.1
PATH=/Library/Frameworks/Python.framework/Versions/Current/bin:${PATH}

# Mecurial
MERCURIAL_DIR=${TOOLS_DIR}/mercurial-0.9.4/${TOOLS_FORMAT}
export PATH=${PATH}:${MERCURIAL_DIR}/bin
export PYTHONPATH=${PYTHONPATH}:${MERCURIAL_DIR}/lib/python${PYTHON_VERSION}/site-packages

# MPICH2
export MPI_DIR=$TOOLS_DIR/mpich2-1.0.5/${TOOLS_FORMAT}
PATH=$PATH:$MPI_DIR/bin
export MANPATH=$MPI_DIR/man:$MANPATH
export MPI_INCDIR=$MPI_DIR/include
export MPI_LIBDIR=$MPI_DIR/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPI_LIBDIR}

# ParMetis
PARMETIS_DIR=$TOOLS_DIR/parmetis-3.1/${TOOLS_FORMAT}
export PARMETIS_INCDIR=$PARMETIS_DIR/include
export PARMETIS_LIBDIR=$PARMETIS_DIR/lib
export PARMETIS_LIBS="-lparmetis -lmetis"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PARMETIS_LIBDIR}

# Cppunit
CPPUNIT_DIR=${TOOLS_DIR}/cppunit-1.12.0/${TOOLS_FORMAT}
PATH=${PATH}:${CPPUNIT_DIR}/bin
export CPPUNIT_LIBDIR=${CPPUNIT_DIR}/lib
export CPPUNIT_INCDIR=${CPPUNIT_DIR}/include
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CPPUNIT_LIBDIR}

# FIAT
FIAT_DIR=${TOOLS_DIR}/fiat-0.3.1/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${FIAT_DIR}/lib/python${PYTHON_VERSION}/site-packages

# PETSc
export PETSC_DIR=${TOOLS_DIR}/petsc-dev
export PETSC_ARCH=osx_${TOOLS_FORMAT}_debug

# Pyrex
PYREX_DIR=${TOOLS_DIR}/pyrex-0.9.5.1/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${PYREX_DIR}/lib/python${PYTHON_VERSION}/site-packages
export PYTHONPATH=${PYTHONPATH}:${PYREX_DIR}/bin
export PATH=${PATH}:${PYREX_DIR}/bin

# Pythia
export PYTHIA_INCDIR=${CIG_INCDIR}/pythia-0.8

# Proj.4
PROJ4_DIR=${TOOLS_DIR}/proj-4.5.0/${TOOLS_FORMAT}
export PROJ4_LIBDIR=${PROJ4_DIR}/lib
export PROJ4_INCDIR=${PROJ4_DIR}/include
PATH=${PATH}:${PROJ4_DIR}/bin
export MANPATH=${MANPATH}:${PROJ4_DIR}/man
