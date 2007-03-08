# -*- autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([testsieve], [0.0], [])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# ----------------------------------------------------------------------
# C/C++/libtool/install
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_INSTALL

# PYTHON
#AM_PATH_PYTHON([2.3])
#CIT_PYTHON_SYSCONFIG

# MPI
AC_LANG(C++)
#CIT_HEADER_MPI
#CIT_CHECK_LIB_MPI

# PETSC
CIT_PATH_PETSC([2.3])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

dnl end of configure.ac
