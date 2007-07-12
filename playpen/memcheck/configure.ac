# -*- autoconf -*-
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([PyLith], [1.0.1], [cig-short@geodynamics.org])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# ----------------------------------------------------------------------
# CUBIT I/O w/netcdf
AC_ARG_ENABLE([cubit],
    [AC_HELP_STRING([--enable-cubit],
        [enable reading/writing Cubit EXODUS files (requires netcdf) @<:@default=no@:>@])],
	[enable_cubit=yes
	CXXFLAGS="-DENABLE_CUBIT $CXXFLAGS"; export CXXFLAGS],
	[enable_cubit=no])
AM_CONDITIONAL([ENABLE_CUBIT], [test "$enable_cubit" = yes])

# PETSc w/TetGen
AC_ARG_ENABLE([tetgen],
    [AC_HELP_STRING([--enable-tetgen],
        [enable generating simple meshes w/TetGen via PETSc @<:@default=no@:>@])],
	[enable_tetgen=yes
	CXXFLAGS="-DENABLE_TETGEN $CXXFLAGS"; export CXXFLAGS],
	[enable_tetgen=no])
AM_CONDITIONAL([ENABLE_TETGEN], [test "$enable_tetgen" = yes])

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
#CIT_PATH_NEMESIS
AM_PATH_PYTHON([2.3])
CIT_PYTHON_SYSCONFIG

# MPI
AC_LANG(C++)
#CIT_PROG_MPICXX
#CIT_HEADER_MPI
#CIT_CHECK_LIB_MPI

# PETSC
CIT_PATH_PETSC([2.3.3])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC
CIT_CHECK_LIB_PETSC_SIEVE

# Check for Python modules and packages.
CIT_PYTHON_EGG_SETUP

# TETGEN
if test "$enable_tetgen" = "yes" ; then
  AC_REQUIRE_CPP
  AC_LANG(C++)
  CPPFLAGS="$PETSC_INCLUDE $PETSC_CXX_INCLUDE $CPPFLAGS"
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

# ENDIANNESS
AC_C_BIGENDIAN

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

dnl end of configure.ac
