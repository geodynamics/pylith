#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $

AC_PREREQ(2.59)
AC_INIT([PyLith3D], [0.8.3], [cig-short@geodynamics.org])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_SRCDIR([pylith3d/pylith3dapp.py.in])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# 'configure' options
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
AC_ARG_ENABLE([embedding],
    [AC_HELP_STRING([--enable-embedding],
        [embed Python with PyLith in a single executable @<:@default=yes@:>@])],
    [want_embedding="$enableval"],
    [want_embedding=yes])
AM_CONDITIONAL([COND_EMBEDDING], [test "$want_embedding" = yes])

# Green's functions option
AC_ARG_ENABLE([greensfns],
    [AC_HELP_STRING([--enable-greensfns],
        [enable Green's functions application (requires Numeric) @<:@default=no@:>@])],
    [enable_greensfns=yes],
    [enable_greensfns=no])
AM_CONDITIONAL([ENABLE_GREENSFNS], [test "$enable_greensfns" = yes])

# Check for Python.
AM_PATH_PYTHON([2.3])
CIT_PYTHON_SYSCONFIG

# Check for Python modules and packages.
CIT_PYTHON_EGG_SETUP

# Check for PETSc.
CIT_PATH_PETSC([2.3.3])
if test -z "$CXX"; then
    AC_PATH_PROG(CXX, $PETSC_CXX, $PETSC_CXX)
    test -z "$CXX" && AC_MSG_FAILURE([C++ compiler $PETSC_CXX (used to configure PETSc) not found in \$PATH])
else
    want_CXX=$CXX; unset CXX
    AC_PATH_PROG(CXX, $want_CXX, $want_CXX)
    test -z "$CXX" && AC_MSG_FAILURE([C++ compiler $want_CXX not found in \$PATH])
fi
if test -z "$FC"; then
    AC_PATH_PROG(FC, $PETSC_FC, $PETSC_FC)
    test -z "$FC" && AC_MSG_FAILURE([Fortran compiler $PETSC_FC (used to configure PETSc) not found in \$PATH])
else
    want_FC=$FC; unset FC
    AC_PATH_PROG(FC, $want_FC, $want_FC)
    test -z "$FC" && AC_MSG_FAILURE([Fortran compiler $want_FC not found in \$PATH])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_FC
F77=$FC
FFLAGS=$FCFLAGS
AC_SUBST(F77)
AC_SUBST(FFLAGS)
AC_PROVIDE([AC_PROG_F77])
AC_FC_WRAPPERS
AC_PROG_LIBTOOL

# Check for MPI.
AC_LANG(C++)
CIT_MPI_CHECK_CXX_LINK(PETSC_INCLUDE, [$PETSC_LIB],
                       _CIT_TRIVIAL_MPI_PROGRAM,
                       [whether we can link a trivial C++ MPI program],
                       [],
                       AC_MSG_FAILURE([cannot link a trivial C++ MPI program]))

# Checks for libraries.
CIT_CHECK_LIB_PETSC
CIT_CHECK_LIB_PETSC_SIEVE

# Checks for header files.
CIT_HEADER_PETSC

# Checks for typedefs, structures, and compiler characteristics.
CIT_FC_STREAM_IO

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 pylith3d/Makefile
                 pylith3d/examples/Makefile
                 pylith3d/libpylith3d/Makefile
                 pylith3d/module/Makefile
		 pylith3d/utils/Makefile])
AC_OUTPUT

dnl end of configure.ac
