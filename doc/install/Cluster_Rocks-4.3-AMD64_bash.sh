# Environment variables for Bash Shells associated with Brad Aagaard's
# install of PyLith on a Linux cluster running Rocks 4.3.

# General
export TOOLS_DIR=${HOME}/tools
export TOOLS_FORMAT=gcc-3.4_64
export PYTHON_VERSION=2.4

# Location for installation of CIG software
export CIG_DIR=${HOME}/tools/cig/${TOOLS_FORMAT}
export CIG_INCDIR=${CIG_DIR}/include
export CIG_LIBDIR=${CIG_DIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CIG_LIBDIR}
export PYTHONPATH=${PYTHONPATH}:${CIG_DIR}/lib64/python${PYTHON_VERSION}/site-packages:${CIG_DIR}/lib/python${PYTHON_VERSION}/site-packages
PATH=${CIG_DIR}/bin:$PATH

# PVFS2
export PVFS2_DIR=/opt/pvfs2
export PVFS2_INCDIR=$PVFS2_DIR/include
export PVFS2_LIBDIR=$PVFS2_DIR/lib

# Python 2.4
PATH=/opt/rocks/bin:${PATH}
export PYTHON_DIR=/opt/rocks/lib/python$PYTHON_VERSION
export PYTHON_LIBDIR=/opt/rocks/lib/python$PYTHON_VERSION
export PYTHON_INCDIR=/opt/rocks/include/python$PYTHON_VERSION

# Numpy
NUMPY_DIR=${TOOLS_DIR}/numpy-1.0.3/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${NUMPY_DIR}/lib/python${PYTHON_VERSION}/site-packages

# Mercurial
MERCURIAL_DIR=${TOOLS_DIR}/mercurial-0.9.3/${TOOLS_FORMAT}
PATH=${PATH}:${MERCURIAL_DIR}/bin
export PYTHONPATH=${PYTHONPATH}:${MERCURIAL_DIR}/lib/python${PYTHON_VERSION}/site-packages

# ACML
ACML_DIR=${TOOLS_DIR}/acml-3.5.0/gnu64
export ACML_LIBDIR=${ACML_DIR}/lib
export ACML_INCDIR=${ACML_DIR}/include
export LD_LIBRARY_PATH=${ACML_LIBDIR}:${LD_LIBRARY_PATH}

# MPICH2
export RSHCOMMAND=ssh
export MPI_DIR=${TOOLS_DIR}/mpich2-1.0.4p1/${TOOLS_FORMAT}
PATH=$MPI_DIR/bin:$PATH
export MANPATH=$MPI_DIR/man:$MANPATH
export MPI_INCDIR=$MPI_DIR/include
export MPI_LIBDIR=$MPI_DIR/lib
export LD_LIBRARY_PATH=${MPI_LIBDIR}:${LD_LIBRARY_PATH}
export MPI_VERSION=2

# Parmetis
PARMETIS_DIR=$TOOLS_DIR/parmetis-3.1/${TOOLS_FORMAT}
export PARMETIS_INCDIR=$PARMETIS_DIR/include
export PARMETIS_LIBDIR=$PARMETIS_DIR/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PARMETIS_LIBDIR}

# Cppunit
CPPUNIT_DIR=${TOOLS_DIR}/cppunit-1.10.2/${TOOLS_FORMAT}
PATH=${PATH}:${CPPUNIT_DIR}/bin
export CPPUNIT_LIBDIR=${CPPUNIT_DIR}/lib
export CPPUNIT_INCDIR=${CPPUNIT_DIR}/include
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CPPUNIT_LIBDIR}

# FIAT
FIAT_DIR=${TOOLS_DIR}/fiat-0.3.3/${TOOLS_FORMAT}
export PYTHONPATH=${PYTHONPATH}:${FIAT_DIR}/lib/python${PYTHON_VERSION}/site-packages

# PETSc
export PETSC_DIR=${TOOLS_DIR}/petsc-dev
export PETSC_ARCH=linux_${TOOLS_FORMAT}_opt

# Pyrex
PYREX_DIR=${TOOLS_DIR}/pyrex-0.9.5.1a/${TOOLS_FORMAT}
PATH=${PATH}:${PYREX_DIR}/bin
export PYTHONPATH=${PYTHONPATH}:${PYREX_DIR}/lib/python${PYTHON_VERSION}/site-packages

# Pythia
export PYTHIA_INCDIR=${CIG_INCDIR}/pythia-0.8

# Proj.4
PROJ4_DIR=${TOOLS_DIR}/proj-4.5.0/${TOOLS_FORMAT}
PATH=${PATH}:${PROJ4_DIR}/bin
export PROJ4_INCDIR=${PROJ4_DIR}/include
export PROJ4_LIBDIR=${PROJ4_DIR}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROJ4_LIBDIR}
export MANPATH=${MANPATH}:${PROJ4_DIR}/man

# NETCDF
NETCDF_DIR=${TOOLS_DIR}/netcdf-3.6.2/${TOOLS_FORMAT}
export NETCDF_INCDIR=${NETCDF_DIR}/include
export NETCDF_LIBDIR=${NETCDF_DIR}/lib
