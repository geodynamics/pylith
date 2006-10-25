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
AC_INIT([PyLith], [1.0.0], [cig-short@geodynamics.org])
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

# TESTING w/cppunit
AC_ARG_ENABLE([testing],
	[  --enable-testing        Enable unit testing with cppunit (requires cppunit) [[default=no]]],
	[enable_testing=yes],
	[enable_testing=no])
AM_CONDITIONAL([ENABLE_TESTING], [test "$enable_testing" = yes])

# DOCUMENTATION w/doxygen
AC_ARG_ENABLE([documentation],
	[  --enable-documentation  Enable buuilding html documentation (requires doxygen) [[default=no]]],
	[enable_documentation=yes],
	[enable_documentation=no])
AM_CONDITIONAL([ENABLE_DOCUMENTATION], [test "$enable_documentation" = yes])

# ----------------------------------------------------------------------
# C/C++/libtool/install
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL

# PYTHON
AM_PATH_PYTHON([2.3])
CIT_PYTHON_SYSCONFIG

# MPI
AC_LANG(C++)
#CIT_PROG_MPICXX
#CIT_HEADER_MPI
#CIT_CHECK_LIB_MPI

# PETSC
CIT_PATH_PETSC([2.3])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

# Checks for Python modules and packages.
CIT_PATH_PYTHIA([0.8], [mpi])
PYTHONPATH="${pyexecdir}${PYTHONPATH:+:${PYTHONPATH}}"; export PYTHONPATH
AC_SUBST([PYTHONPATH])

# PYREX/PYREXEMBED
AC_CHECK_PROG([have_pyrex],
  [pyrexc],
  [yes],
  [no])
if test "$have_pyrex" = "no"; then
  AC_MSG_ERROR([Pyrex not found in current path.])
fi
AC_CHECK_PROG([have_pyrexembed],
  [pyrexembed],
  [yes],
  [no])
if test "$have_pyrexembed" = "no"; then
  AC_MSG_ERROR([Pyrexembed not found in current path.])
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

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
		pylith/Makefile
                libsrc/Makefile
                libsrc/feassemble/Makefile
                libsrc/meshio/Makefile
		libsrc/utils/Makefile
                modulesrc/Makefile
                modulesrc/feassemble/Makefile
                modulesrc/meshio/Makefile
		applications/Makefile
		unittests/Makefile
		unittests/libtests/Makefile
		unittests/libtests/feassemble/Makefile
		unittests/pytests/Makefile
		unittests/pytests/feassemble/Makefile
                doc/Makefile])

#                 libsrc/materials/Makefile

AC_OUTPUT

dnl end of configure.ac
