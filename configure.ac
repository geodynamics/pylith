# -*- autoconf -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([PyLith], [2.1.0], [cig-short@geodynamics.org])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])

# ----------------------------------------------------------------------
# PYTHON
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])

# SWIG
AC_ARG_ENABLE([swig],
    [AC_HELP_STRING([--enable-swig],
        [enable generating modules with SWIG) @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then enable_swig=yes; else enable_swig=no; fi],
	[enable_swig=no])
AM_CONDITIONAL([ENABLE_SWIG], [test "$enable_swig" = yes])

# TESTING w/cppunit
AC_ARG_ENABLE([testing],
    [AC_HELP_STRING([--enable-testing],
        [enable unit testing with cppunit (requires cppunit) @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then enable_testing=yes; else enable_testing=no; fi],
	[enable_testing=no])
AM_CONDITIONAL([ENABLE_TESTING], [test "$enable_testing" = yes])
# Enable full-scale testing except on cygwin (need to rebase)
if test "$enable_testing" = yes ; then
  enable_full_testing=yes
  AC_CANONICAL_HOST
  case $host_os in
    cygwin*)
      enable_full_testing=no ;;
  esac
fi
AM_CONDITIONAL([ENABLE_FULL_TESTING], [test "$enable_full_testing" = yes])

# CUBIT I/O w/netcdf
AC_ARG_ENABLE([cubit],
    [AC_HELP_STRING([--enable-cubit],
        [enable reading/writing Cubit EXODUS files (requires netcdf) @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then 
	  enable_cubit=yes
	  CPPFLAGS="-DENABLE_CUBIT $CPPFLAGS"; export CPPFLAGS;
	  PYLITH_SWIG_CPPFLAGS="-DENABLE_CUBIT $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS
	else enable_cubit=no; fi],
	[enable_cubit=no])
AM_CONDITIONAL([ENABLE_CUBIT], [test "$enable_cubit" = yes])
AC_SUBST(PYLITH_SWIG_CPPFLAGS)

# PETSc w/TetGen
AC_ARG_ENABLE([tetgen],
    [AC_HELP_STRING([--enable-tetgen],
        [enable generating simple meshes w/TetGen via PETSc @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then 
	  enable_tetgen=yes
	  CPPFLAGS="-DENABLE_TETGEN $CPPFLAGS"; export CPPFLAGS
        else enable_tetgen=no; fi],
	[enable_tetgen=no])
AM_CONDITIONAL([ENABLE_TETGEN], [test "$enable_tetgen" = yes])

# PETSc w/HDF5
AC_ARG_ENABLE([hdf5],
    [AC_HELP_STRING([--enable-hdf5],
        [enable output to HDF5 files via PETSc @<:@default=no@:>@])],
	[if test "$enableval" = yes; 
	  then enable_hdf5=yes
          else enable_hdf5=no; fi],
	[enable_hdf5=yes])
if test "$enable_hdf5" = "yes"; then
  CPPFLAGS="-DENABLE_HDF5 $CPPFLAGS"; export CPPFLAGS
  PYLITH_SWIG_CPPFLAGS="-DENABLE_HDF5 $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS
fi
AM_CONDITIONAL([ENABLE_HDF5], [test "$enable_hdf5" = yes])

# CUDA
AC_ARG_ENABLE([cuda],
    [AC_HELP_STRING([--enable-cuda],
        [enable use of CUDA @<:@default=no@:>@])],
	[if test "$enable_cuda" = yes; then
          enable_cuda=yes
          CPPFLAGS="-DENABLE_CUDA $CPPFLAGS"; export CPPFLAGS
	  PYLITH_SWIG_CPPFLAGS="-DENABLE_CUDA $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS
          else enable_cuda=no; fi],
	[enable_cuda=no])
AM_CONDITIONAL([ENABLE_CUDA], [test "$enable_cuda" = yes])

# MEMORY LOGGING
AC_ARG_ENABLE([memory_logging],
    [AC_HELP_STRING([--enable-memory-logging],
        [enable memory logging @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then 
	  enable_memory_logging=yes
	  CPPFLAGS="-DALE_MEM_LOGGING $CPPFLAGS"; export CPPFLAGS
         else enable_memory_logging=no; fi],
	[enable_memory_logging=no])
AM_CONDITIONAL([ENABLE_MEMORY_LOGGING], [test "$enable_memory_logging" = yes])

# DOCUMENTATION w/doxygen
AC_ARG_ENABLE([documentation],
    [AC_HELP_STRING([--enable-documentation],
        [enable building HTML documentation (requires doxygen) @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then enable_documentation=yes; else enable_documentation=no; fi],
	[enable_documentation=no])
AM_CONDITIONAL([ENABLE_DOCUMENTATION], [test "$enable_documentation" = yes])

# ----------------------------------------------------------------------
# C/C++/libtool/install
AC_PROG_CXX
AC_PROG_CC
AC_DISABLE_STATIC

AC_PROG_LIBTOOL
if test "$allow_undefined_flag" = unsupported; then
    # See issue119.
    AM_LDFLAGS="-no-undefined $AM_LDFLAGS"
fi
AM_CONDITIONAL([NO_UNDEFINED], [test "$allow_undefined_flag" = unsupported])
AC_SUBST(AM_LDFLAGS)

AC_PROG_INSTALL

# PYTHON
CIT_PATH_NEMESIS
AM_PATH_PYTHON([2.7])
CIT_PYTHON_SYSCONFIG
CIT_CHECK_PYTHON_HEADER
CIT_CHECK_PYTHON_SHARED

# MPI
AC_LANG(C++)
#CIT_PROG_MPICXX
#CIT_HEADER_MPI
#CIT_CHECK_LIB_MPI

# PETSC
AC_LANG(C)
CIT_PATH_PETSC([3.6.0])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

# Large file support
AC_SYS_LARGEFILE

# NUMPY
CIT_NUMPY_PYTHON_MODULE
CIT_NUMPY_INCDIR

# FIAT
CIT_PYTHON_MODULE([FIAT],[0.9.9])

# SWIG
if test "$enable_swig" = "yes" ; then
  AC_PROG_SWIG(1.3.33)
fi

# CPPUNIT
if test "$enable_testing" = "yes" ; then
  CIT_CPPUNIT_HEADER
  CIT_CPPUNIT_LIB
fi

# CUBIT (netcdf)
if test "$enable_cubit" = "yes" ; then
  CIT_NETCDF_HEADER
  CIT_NETCDF_LIB
fi

# TETGEN
if test "$enable_tetgen" = "yes" ; then
  AC_REQUIRE_CPP
  AC_LANG(C++)
  tetgen_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PETSC_CC_INCLUDES $PETSC_CXX_INCLUDE $CPPFLAGS"
  AC_CHECK_HEADER([tetgen.h], [], [
    AC_MSG_ERROR([tetgen header not found; check PETSc installation.])
  ])
  AC_MSG_CHECKING([for tetgenbehavior in -ltetgen])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <tetgen.h>]],
	             [[tetgenbehavior  behavior;]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([tetgen library not found; check PETSc installation.])
  ])
  CPPFLAGS=$tetgen_save_CPPFLAGS
fi

# HDF5
if test "$enable_hdf5" = "yes" ; then
  CIT_HDF5_HEADER
  CIT_HDF5_LIB
  CIT_HDF5_LIB_PARALLEL
fi

# Full-scale testing with h5py
if test "$enable_full_testing" = yes ; then
  AM_PATH_PYTHON
   AC_MSG_CHECKING(for python module h5py for full scale testing)
  $PYTHON -c "import h5py" 2>/dev/null   
  if test $? == 0; then
    AC_MSG_RESULT([found h5py, enabling full scale testing])
    enable_full_testing=yes
  else
    AC_MSG_RESULT([not found, disabling full scale testing])
    enable_full_testing=no
  fi
fi
AM_CONDITIONAL([ENABLE_FULL_TESTING], [test "$enable_full_testing" = yes])

# CUDA
if test "$enable_cuda" = "yes" ; then
  CIT_CUDA_INCDIR
  CIT_CUDA_LIB
  CIT_CUDA_COMPILER
fi

# PROJ4
CIT_PROJ4_HEADER
CIT_PROJ4_LIB

# PYTHIA
CIT_PYTHON_MODULE([pyre],[0.8.1])
CIT_PYTHON_MODULE([journal],[0.8.1])
eval p=`$PYTHON -c "import pyre; print pyre.__""file__.replace('pyre/__init__.pyc','')"`
PYTHON_EGG_CPPFLAGS="-I$p"
AC_SUBST(PYTHON_EGG_CPPFLAGS)

# SPATIALDATA
CIT_SPATIALDATA_HEADER
CIT_SPATIALDATA_LIB
CIT_PYTHON_MODULE([spatialdata],[1.9.3])

# DOXYGEN
if test "$enable_documentation" = "yes" ; then
  AC_CHECK_PROG([have_doxygen],
    [doxygen],
    [yes],
    [no])
  if test "$have_doxygen" = "no"; then
    AC_MSG_ERROR([Doxygen not found in current path.])
  fi
fi

# Check whether compiler sets variable with function names.
CIT_FUNCTIONSTRING

# ENDIANNESS
AC_C_BIGENDIAN

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
		pylith/Makefile
		libsrc/Makefile
		libsrc/pylith/Makefile
		libsrc/pylith/bc/Makefile
                libsrc/pylith/feassemble/Makefile
		libsrc/pylith/faults/Makefile
		libsrc/pylith/friction/Makefile
                libsrc/pylith/materials/Makefile
                libsrc/pylith/meshio/Makefile
                libsrc/pylith/problems/Makefile
		libsrc/pylith/topology/Makefile
		libsrc/pylith/utils/Makefile
                modulesrc/Makefile
                modulesrc/include/Makefile
                modulesrc/bc/Makefile
                modulesrc/faults/Makefile
                modulesrc/feassemble/Makefile
		modulesrc/friction/Makefile
                modulesrc/materials/Makefile
                modulesrc/meshio/Makefile
                modulesrc/mpi/Makefile
                modulesrc/problems/Makefile
		modulesrc/topology/Makefile
		modulesrc/utils/Makefile
		applications/Makefile
		applications/utilities/Makefile
		unittests/Makefile
		unittests/libtests/Makefile
		unittests/libtests/bc/Makefile
		unittests/libtests/bc/data/Makefile
		unittests/libtests/faults/Makefile
		unittests/libtests/faults/data/Makefile
		unittests/libtests/feassemble/Makefile
		unittests/libtests/feassemble/data/Makefile
		unittests/libtests/friction/Makefile
		unittests/libtests/friction/data/Makefile
		unittests/libtests/materials/Makefile
		unittests/libtests/materials/data/Makefile
		unittests/libtests/meshio/Makefile
		unittests/libtests/meshio/data/Makefile
		unittests/libtests/topology/Makefile
		unittests/libtests/topology/data/Makefile
		unittests/libtests/utils/Makefile
		unittests/pytests/Makefile
		unittests/pytests/bc/Makefile
		unittests/pytests/bc/data/Makefile
		unittests/pytests/faults/Makefile
		unittests/pytests/faults/data/Makefile
		unittests/pytests/feassemble/Makefile
		unittests/pytests/feassemble/data/Makefile
		unittests/pytests/friction/Makefile
		unittests/pytests/friction/data/Makefile
		unittests/pytests/materials/Makefile
		unittests/pytests/materials/data/Makefile
		unittests/pytests/meshio/Makefile
		unittests/pytests/meshio/data/Makefile
		unittests/pytests/mpi/Makefile
		unittests/pytests/problems/Makefile
		unittests/pytests/problems/data/Makefile
		unittests/pytests/topology/Makefile
		unittests/pytests/topology/data/Makefile
		unittests/pytests/utils/Makefile
		tests_auto/Makefile
		tests_auto/2d/Makefile
		tests_auto/2d/tri3/Makefile
		tests_auto/2d/quad4/Makefile
		tests_auto/2d/slipdir/Makefile
		tests_auto/3d/Makefile
		tests_auto/3d/tet4/Makefile
		tests_auto/3d/hex8/Makefile
		tests_auto/petsc/Makefile
		tests_auto/eqinfo/Makefile
		tests/Makefile
		tests/2d/Makefile
		tests/2d/maxwell/Makefile
		tests/2d/powerlaw/Makefile
		tests/2d/frictionslide/Makefile
		tests/2d/plasticity/Makefile
		tests/2d/plasticity/initialstress/Makefile
		tests/3d/Makefile
		tests/3d/matprops/Makefile
		tests/3d/slipdir/Makefile
		tests/3d/cyclicfriction/Makefile
		tests/3d/cyclicfriction/output/Makefile
		tests/3d/plasticity/Makefile
		tests/3d/plasticity/dynamic/Makefile
		tests/3d/plasticity/initialstress/Makefile
                doc/Makefile
		doc/developer/Makefile
		doc/install/Makefile
		doc/refguide/Makefile
		doc/releasenotes/Makefile
		doc/userguide/Makefile
		doc/presentations/Makefile
		examples/Makefile
		examples/3d/Makefile
		examples/3d/tet4/Makefile
		examples/3d/tet4/output/Makefile
		examples/3d/hex8/Makefile
		examples/3d/hex8/output/Makefile
		examples/2d/Makefile
		examples/2d/subduction/Makefile
		examples/2d/subduction/output/Makefile
		examples/2d/gravity/Makefile
		examples/2d/greensfns/Makefile
		examples/2d/greensfns/strikeslip/Makefile
		examples/2d/greensfns/strikeslip/output/Makefile
		examples/2d/greensfns/reverse/Makefile
		examples/2d/greensfns/reverse/output/Makefile
		examples/bar_shearwave/Makefile
		examples/bar_shearwave/hex8/Makefile
		examples/bar_shearwave/hex8/output/Makefile
		examples/bar_shearwave/quad4/Makefile
		examples/bar_shearwave/quad4/output/Makefile
		examples/bar_shearwave/tet4/Makefile
		examples/bar_shearwave/tet4/output/Makefile
		examples/bar_shearwave/tri3/Makefile
		examples/bar_shearwave/tri3/output/Makefile
		examples/twocells/Makefile
		examples/twocells/twohex8/Makefile
		examples/twocells/twoquad4/Makefile
		examples/twocells/twotet4/Makefile
		examples/twocells/twotet4-geoproj/Makefile
		examples/twocells/twotri3/Makefile
		examples/meshing/Makefile
		examples/meshing/surface_nurbs/Makefile
		examples/meshing/surface_nurbs/dem/Makefile
		examples/meshing/surface_nurbs/dem/ulines/Makefile
		examples/meshing/surface_nurbs/dem/vlines/Makefile
		examples/meshing/surface_nurbs/triangles/Makefile
		examples/meshing/surface_nurbs/merge_surfs/Makefile
		examples/meshing/surface_nurbs/subduction/Makefile
		examples/meshing/cubit_cellsize/Makefile
		examples/debugging/Makefile
		templates/Makefile
		share/Makefile
		])

AC_OUTPUT

dnl end of configure.ac
