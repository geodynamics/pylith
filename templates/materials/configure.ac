# -*- autoconf -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

AC_PREREQ(2.59)
AC_INIT([materialscontrib], [0.0.1], [cig-short@geodynamics.org])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_MACRO_DIR([./m4])
AM_INIT_AUTOMAKE([foreign])

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
AM_PATH_PYTHON([2.7])
CIT_PYTHON_SYSCONFIG
CIT_CHECK_PYTHON_HEADER
CIT_CHECK_PYTHON_SHARED

# SWIG
CIT_NUMPY_PYTHON_MODULE
CIT_NUMPY_INCDIR
AC_PROG_SWIG([2.0.4],[3.0.2])

# MPI
AC_LANG(C++)
CIT_HEADER_MPI

# PETSC (dictates precision)
AC_LANG(C++)
CIT_PATH_PETSC([3.8.2])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

# PYLITH
pylith_save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$PETSC_CC_INCLUDES $CPPFLAGS"
AC_LANG(C++)
AC_CHECK_HEADER([pylith/materials/ElasticMaterial.hh], [], [
  AC_MSG_ERROR([PyLith ElasticMaterial header not found; try CPPFLAGS="-I<PyLith include dir>"])
])
AC_MSG_CHECKING([for materials::ElasticIsotropic3D in -lpylith])
AC_REQUIRE_CPP
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([[#include <pylith/materials/ElasticIsotropic3D.hh>]],
                   [[pylith::materials::ElasticIsotropic3D material;]])],
  [AC_MSG_RESULT(yes)],
  [AC_MSG_RESULT(no)
   AC_MSG_ERROR([PyLith library not found; try LDFLAGS="-L<PyLith lib dir>"])
])
if test "x$PYLITH_SWIG_DIR" != "x" ; then
  AC_CHECK_FILE([$PYLITH_SWIG_DIR/materials/ElasticMaterial.i], [], [
    AC_MSG_ERROR([PyLith ElasticMaterial.i SWIG interface file not found])])
  AC_SUBST([SWIG_FLAGS], ["-I$PYLITH_SWIG_DIR $SWIG_FLAGS"])
else
  AC_CHECK_FILE([materials/ElasticMaterial.i], [], [
    AC_MSG_ERROR([PyLith ElasticMaterial.i SWIG interface file not found; Try setting PYLITH_SWIG_DIR=<directory containing materials/ElasticMaterial.i>])])
fi
CPPFLAGS=$pylith_save_CPPFLAGS


# ENDIANNESS
AC_C_BIGENDIAN

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
	         tests/Makefile])

AC_OUTPUT


# End of file
