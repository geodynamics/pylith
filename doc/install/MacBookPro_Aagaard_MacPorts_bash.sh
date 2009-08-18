# Environment variables for Bash Shells associated with Brad Aagaard's
# install of PyLith on a MacBookPro running OS X 10.5

# General
export TOOLS_DIR=/tools
export TOOLS_FORMAT=gcc-4.0
export PYTHON_VERSION=2.5

# Location for installation of CIG software
export CIG_DIR=$HOME/tools/cig/${TOOLS_FORMAT}
export CIG_INCDIR=${CIG_DIR}/include
export CIG_LIBDIR=${CIG_DIR}/lib
PATH=$PATH:${CIG_DIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CIG_LIBDIR}
export PYTHONPATH=${PYTHONPATH}:${CIG_LIBDIR}/python${PYTHON_VERSION}/site-packages

# MPICH2
export MPI_DIR=$TOOLS_DIR/mpich2-1.1.1/${TOOLS_FORMAT}
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

# FIAT
FIAT_DIR=${TOOLS_DIR}/fiat-0.3.4/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${FIAT_DIR}/lib/python${PYTHON_VERSION}/site-packages

# PETSc
export PETSC_DIR=${TOOLS_DIR}/petsc-dev
export PETSC_ARCH=osx_${TOOLS_FORMAT}_debug

# Pythia
export PYTHIA_INCDIR=${CIG_INCDIR}/pythia-0.8

# Proj.4
PROJ4_DIR=${TOOLS_DIR}/proj-4.6.1/${TOOLS_FORMAT}
export PROJ4_LIBDIR=${PROJ4_DIR}/lib
export PROJ4_INCDIR=${PROJ4_DIR}/include
PATH=${PATH}:${PROJ4_DIR}/bin
export MANPATH=${MANPATH}:${PROJ4_DIR}/man
