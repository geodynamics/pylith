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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([PyLith], [1.5.1], [cig-short@geodynamics.org])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# ----------------------------------------------------------------------
# PYTHON
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
AC_ARG_WITH([embedding],
    [AC_HELP_STRING([--with-embedding],
        [embed Python with PyLith in a single executable @<:@default=yes@:>@])],
    [want_embedding="$withval"],
    [want_embedding=no])
AM_CONDITIONAL([COND_EMBEDDING], [test "$want_embedding" = yes])

# SWIG
AC_ARG_ENABLE([swig],
    [AC_HELP_STRING([--enable-swig],
        [enable generating modules with SWIG) @<:@default=no@:>@])],
	[enable_swig=yes],
	[enable_swig=no])
AM_CONDITIONAL([ENABLE_SWIG], [test "$enable_swig" = yes])

# TESTING w/cppunit
AC_ARG_ENABLE([testing],
    [AC_HELP_STRING([--enable-testing],
        [enable unit testing with cppunit (requires cppunit) @<:@default=no@:>@])],
	[enable_testing=yes],
	[enable_testing=no])
AM_CONDITIONAL([ENABLE_TESTING], [test "$enable_testing" = yes])

# CUBIT I/O w/netcdf
AC_ARG_ENABLE([cubit],
    [AC_HELP_STRING([--enable-cubit],
        [enable reading/writing Cubit EXODUS files (requires netcdf) @<:@default=no@:>@])],
	[enable_cubit=yes
	CPPFLAGS="-DENABLE_CUBIT $CPPFLAGS"; export CPPFLAGS;
	PYLITH_SWIG_CPPFLAGS="-DENABLE_CUBIT $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS],
	[enable_cubit=no])
AM_CONDITIONAL([ENABLE_CUBIT], [test "$enable_cubit" = yes])
AC_SUBST(PYLITH_SWIG_CPPFLAGS)

# PETSc w/TetGen
AC_ARG_ENABLE([tetgen],
    [AC_HELP_STRING([--enable-tetgen],
        [enable generating simple meshes w/TetGen via PETSc @<:@default=no@:>@])],
	[enable_tetgen=yes
	CPPFLAGS="-DENABLE_TETGEN $CPPFLAGS"; export CPPFLAGS],
	[enable_tetgen=no])
AM_CONDITIONAL([ENABLE_TETGEN], [test "$enable_tetgen" = yes])

# PETSc w/HDF5
AC_ARG_ENABLE([hdf5],
    [AC_HELP_STRING([--enable-hdf5],
        [enable output to HDF5 files via PETSc @<:@default=no@:>@])],
	[enable_hdf5=yes
	CPPFLAGS="-DENABLE_HDF5 $CPPFLAGS"; export CPPFLAGS
	PYLITH_SWIG_CPPFLAGS="-DENABLE_HDF5 $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS],
	[enable_hdf5=no])
AM_CONDITIONAL([ENABLE_HDF5], [test "$enable_hdf5" = yes])

# CUDA
AC_ARG_ENABLE([cuda],
    [AC_HELP_STRING([--enable-cuda],
        [enable use of CUDA @<:@default=no@:>@])],
	[enable_cuda=yes
	CPPFLAGS="-DENABLE_CUDA $CPPFLAGS"; export CPPFLAGS
	PYLITH_SWIG_CPPFLAGS="-DENABLE_CUDA $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS],
	[enable_cuda=no])
AM_CONDITIONAL([ENABLE_CUDA], [test "$enable_cuda" = yes])

# MEMORY LOGGING
AC_ARG_ENABLE([memory_logging],
    [AC_HELP_STRING([--enable-memory-logging],
        [enable memory logging @<:@default=no@:>@])],
	[enable_memory_logging=yes
	CPPFLAGS="-DALE_MEM_LOGGING $CPPFLAGS"; export CPPFLAGS],
	[enable_memory_logging=no])
AM_CONDITIONAL([ENABLE_MEMORY_LOGGING], [test "$enable_memory_logging" = yes])

# DOCUMENTATION w/doxygen
AC_ARG_ENABLE([documentation],
    [AC_HELP_STRING([--enable-documentation],
        [enable building HTML documentation (requires doxygen) @<:@default=no@:>@])],
	[enable_documentation=yes],
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
AM_PATH_PYTHON([2.3])
CIT_PYTHON_SYSCONFIG

# MPI
AC_LANG(C++)
#CIT_PROG_MPICXX
#CIT_HEADER_MPI
#CIT_CHECK_LIB_MPI

# PETSC
CIT_PATH_PETSC([3.1.0])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC
CIT_CHECK_LIB_PETSC_SIEVE

# Check for Python modules and packages.
CIT_PYTHON_EGG_SETUP

# NUMPY
CIT_NUMPY_PYTHON_MODULE
CIT_NUMPY_INCDIR

# SWIG
if test "$enable_swig" = "yes" ; then
  AC_PROG_SWIG(1.3.33)
fi

# CPPUNIT
if test "$enable_testing" = "yes" ; then
  AC_LANG(C++)
  AC_CHECK_HEADER([cppunit/TestRunner.h], [], [
    AC_MSG_ERROR([CppUnit header not found; try CPPFLAGS="-I<CppUnit include dir>"])
  ])
  AC_MSG_CHECKING([for CppUnit::TestRunner in -lcppunit])
  AC_REQUIRE_CPP
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <cppunit/TestRunner.h>]],
	             [[CppUnit::TestRunner runner;]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([CppUnit library not found; try LDFLAGS="-L<CppUnit lib dir>"])
  ])
fi

# CUBIT (netcdf)
if test "$enable_cubit" = "yes" ; then
  AC_LANG(C++)
  AC_CHECK_HEADER([netcdfcpp.h], [], [
    AC_MSG_ERROR([netcdf C++ header not found; try CPPFLAGS="-I<netcdf include dir>"])
  ])
  AC_MSG_CHECKING([for NcFile in -lnetcdfc++])
  AC_REQUIRE_CPP
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <netcdfcpp.h>]],
	             [[NcFile ncfile("filename");]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([netcdfc++ library not found; try LDFLAGS="-L<netcdf lib dir>"])
  ])
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
  AC_REQUIRE_CPP
  AC_LANG(C++)
  hdf5_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PETSC_CC_INCLUDES $PETSC_CXX_INCLUDE $CPPFLAGS"  
  AC_CHECK_HEADER([hdf5.h], [], [
    AC_MSG_ERROR([HDF5 header not found.])
  ])
  AC_MSG_CHECKING([for H5Fopen in -lhdf5])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <hdf5.h>]],
	             [[H5Fopen("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT);]])],
    [AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)
     AC_MSG_ERROR([hdf5 library not found.])
  ])
  CPPFLAGS=$hdf5_save_CPPFLAGS
fi

# CUDA
if test "$enable_cuda" = "yes" ; then
  AC_REQUIRE_CPP
  AC_LANG(C++)
  cuda_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PETSC_CC_INCLUDES $PETSC_CXX_INCLUDE $CPPFLAGS"
  AC_CHECK_HEADER([cuda_runtime.h], [], [
    AC_MSG_ERROR([CUDA runtime header not found.])
  ])
  #AC_MSG_CHECKING([for ?? in -lcuda])
  #AC_COMPILE_IFELSE(
  #  [AC_LANG_PROGRAM([[#include <cuda_runtime.h>]],
  #	             [[H5Fopen("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT);]])],
  #  [AC_MSG_RESULT(yes)],
  #  [AC_MSG_RESULT(no)
  #   AC_MSG_ERROR([cuda library not found.])
  #])
  CPPFLAGS=$cuda_save_CPPFLAGS
fi

# PROJ4
AC_CHECK_LIB(proj, pj_init_plus, [
  AC_CHECK_HEADER([proj_api.h], [], [
    AC_MSG_ERROR([Proj4 header not found; try CPPFLAGS="-I<Proj4 include dir>"])
  ])
],[
  AC_MSG_ERROR([Proj4 library not found; try LDFLAGS="-L<Proj4 lib dir>"])
])

# SPATIALDATA
AC_LANG(C++)
AC_CHECK_HEADER([spatialdata/spatialdb/SpatialDB.hh], [], [
  AC_MSG_ERROR([SpatialDB header not found; try CPPFLAGS="-I<Spatialdata include dir>"])
])
AC_MSG_CHECKING([for spatialdb::SimpleDB in -lspatialdata])
AC_REQUIRE_CPP
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <spatialdata/spatialdb/SpatialDB.hh>]
                    [#include <spatialdata/spatialdb/SimpleDB.hh>]],
                   [[spatialdata::spatialdb::SimpleDB db;]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_MSG_ERROR([Spatialdata library not found; try LDFLAGS="-L<Spatialdata lib dir>"])
])

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
		libsrc/bc/Makefile
                libsrc/feassemble/Makefile
		libsrc/faults/Makefile
		libsrc/friction/Makefile
                libsrc/materials/Makefile
                libsrc/meshio/Makefile
                libsrc/problems/Makefile
		libsrc/topology/Makefile
		libsrc/utils/Makefile
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
		tests_auto/1d/Makefile
		tests_auto/1d/line2/Makefile
		tests_auto/1d/line3/Makefile
		tests_auto/2d/Makefile
		tests_auto/2d/tri3/Makefile
		tests_auto/2d/quad4/Makefile
		tests_auto/petsc/Makefile
		tests/Makefile
		tests/2d/Makefile
		tests/2d/maxwell/Makefile
		tests/2d/slipdir/Makefile
		tests/3d/Makefile
		tests/3d/matprops/Makefile
		tests/3d/slipdir/Makefile
		tests/3d/plasticity/Makefile
		tests/3d/plasticity/threehex27/Makefile
		tests/3d/plasticity/threehex27/config/Makefile
		tests/3d/plasticity/threehex27/mesh/Makefile
		tests/3d/plasticity/threehex27/results/Makefile
		tests/3d/plasticity/threehex27/spatialdb/Makefile
		tests/3d/plasticity/threehex8/Makefile
		tests/3d/plasticity/threehex8/config/Makefile
		tests/3d/plasticity/threehex8/mesh/Makefile
		tests/3d/plasticity/threehex8/results/Makefile
		tests/3d/plasticity/threehex8/spatialdb/Makefile
		tests/topology/Makefile
                doc/Makefile
		doc/developer/Makefile
		doc/install/Makefile
		doc/refguide/Makefile
		doc/releasenotes/Makefile
		doc/userguide/Makefile
		doc/presentations/Makefile
		examples/Makefile
		examples/3d/Makefile
		examples/bar_shearwave/Makefile
		examples/twocells/Makefile
		examples/3d/tet4/Makefile
		examples/3d/hex8/Makefile
		examples/3d/hex8/output/Makefile
		examples/bar_shearwave/hex8/Makefile
		examples/bar_shearwave/hex8/output/Makefile
		examples/bar_shearwave/quad4/Makefile
		examples/bar_shearwave/quad4/output/Makefile
		examples/bar_shearwave/tet4/Makefile
		examples/bar_shearwave/tet4/output/Makefile
		examples/bar_shearwave/tri3/Makefile
		examples/bar_shearwave/tri3/output/Makefile
		examples/greensfns/Makefile
		examples/greensfns/hex8/Makefile
		examples/greensfns/hex8/gfimpulses/Makefile
		examples/greensfns/hex8/gfresponses/Makefile
		examples/greensfns/hex8/gfspatialdb/Makefile
		examples/twocells/twohex8/Makefile
		examples/twocells/twohex27-cubit/Makefile
		examples/twocells/twoquad4/Makefile
		examples/twocells/twotet4/Makefile
		examples/twocells/twotet4-geoproj/Makefile
		examples/twocells/twotri3/Makefile
		examples/meshing/Makefile
		examples/meshing/surface_nurbs/Makefile
		examples/meshing/surface_nurbs/contours/Makefile
		examples/meshing/surface_nurbs/contours/journal/Makefile
		examples/meshing/surface_nurbs/dem/Makefile
		examples/meshing/surface_nurbs/dem/ulines/Makefile
		examples/meshing/surface_nurbs/dem/vlines/Makefile
		examples/meshing/surface_nurbs/triangles/Makefile
		templates/Makefile
		share/Makefile
		])

AC_OUTPUT

dnl end of configure.ac
