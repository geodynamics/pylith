# Environment variables for Bash Shells associated with Brad Aagaard's
# install of PyLith on a AMD64 computer running Fedora 8.

# General
export TOOLS_DIR=/tools/common
export TOOLS_FORMAT=gcc-4.1.2_64
export PYTHON_VERSION=2.5

# Location for installation of CIG software
export CIG_DIR=${HOME}/tools/cig/${TOOLS_FORMAT}
export CIG_INCDIR=${CIG_DIR}/include
export CIG_LIBDIR=${CIG_DIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CIG_LIBDIR}
export PYTHONPATH=${PYTHONPATH}:${CIG_DIR}/lib64/python${PYTHON_VERSION}/site-packages:${CIG_DIR}/lib/python${PYTHON_VERSION}/site-packages
PATH=${CIG_DIR}/bin:$PATH

# ACML
ACML_DIR=${TOOLS_DIR}/acml-4.0.1/gfortran64
export ACML_INCDIR=${ACML_DIR}/include
export ACML_LIBDIR=${ACML_DIR}/lib
export ACML_LIBS="-lacml_mv -lacml -lgfortran -lm"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ACML_LIBDIR}

# MPICH2
export RSHCOMMAND=ssh
export MPI_DIR=$TOOLS_DIR/mpich2-1.0.6p1/${TOOLS_FORMAT}
PATH=$MPI_DIR/bin:$PATH
export MANPATH=$MPI_DIR/man:$MANPATH
export MPI_INCDIR=$MPI_DIR/include
export MPI_LIBDIR=$MPI_DIR/lib
export LD_LIBRARY_PATH=${MPI_LIBDIR}:${LD_LIBRARY_PATH}

# PARMETIS
PARMETIS_DIR=$TOOLS_DIR/parmetis-3.1/${TOOLS_FORMAT}
export PARMETIS_INCDIR=$PARMETIS_DIR/include
export PARMETIS_LIBDIR=$PARMETIS_DIR/lib
export PARMETIS_LIBS="-lparmetis -lmetis"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PARMETIS_LIBDIR}

# FIAT
FIAT_DIR=${TOOLS_DIR}/fiat-0.3.1/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${FIAT_DIR}/lib/python${PYTHON_VERSION}/site-packages

# PETSc
export PETSC_DIR=${TOOLS_DIR}/petsc-dev
export PETSC_ARCH=linux_${TOOLS_FORMAT}_debug

# Pythia
export PYTHIA_INCDIR=${CIG_INCDIR}/pythia-0.8

# NETCDF (installed as package)
export NETCDF_INCDIR=/usr/include/netcdf-3
export NETCDF_LIBDIR=/usr/lib64/netcdf-3
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETCDF_LIBDIR}
