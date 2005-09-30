#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $

AC_PREREQ(2.59)
AC_INIT([Lithomop3d], [1.0.0], [willic3@rpi.edu])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_SRCDIR([lithomop3d/applications/lithomop3dapp.py])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# 'configure' options
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])

# Checks for programs.
AM_PATH_PYTHON([2.3])
CIT_PYTHON_INCDIR
CIT_PATH_PETSC([2.3])
CIT_PROG_MPICXX
# this is not portable
CXX=$MPICXX
if test -z "$FC"; then
    AC_CHECK_PROG(FC, $PETSC_FC, $PETSC_FC)
    test -z "$FC" && AC_MSG_FAILURE([Fortran compiler $PETSC_FC (used to configure PETSc) not found in \$PATH])
else
    want_FC=$FC; unset FC
    AC_CHECK_PROG(FC, $want_FC, $want_FC)
    test -z "$FC" && AC_MSG_FAILURE([Fortran compiler $want_FC not found in \$PATH])
fi
AC_PROG_FC
AC_FC_WRAPPERS
AC_PROG_LIBTOOL

# Checks for Python modules and packages.
CIT_PATH_PYTHIA([0.8], [mpi])

# Checks for libraries.
CIT_CHECK_LIB_MPI

# Checks for header files.
CIT_HEADER_MPI

# Checks for typedefs, structures, and compiler characteristics.
AC_LANG(Fortran)
AC_MSG_CHECKING([whether $FC supports OPEN control item 'position'])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(kw,file="foo",status="old",position="append")]])
], [
    AC_MSG_RESULT(yes)
], [
    AC_MSG_RESULT(no)
    AC_MSG_CHECKING([whether $FC supports OPEN control item 'access'])
    AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [[      open(kw,file="foo",status="old",access="append")]])
    ], [
        AC_MSG_RESULT(yes)
        FCFLAGS="-DFORTRAN_TYPE_F77 $FCFLAGS"; export FCFLAGS
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_FAILURE([can not determine Fortran type])
    ])
])
AC_MSG_CHECKING([whether $FC supports stream I/O])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[
       integer foo,kucd
       open(kucd,file="foo",status="new",access="stream",
     &   form="unformatted")
       write(kucd,pos=ibyte) foo
]])], [AC_MSG_RESULT(yes)],
      [AC_MSG_RESULT(no)
       FCFLAGS="-DFORTRAN_NAME_GFORTRAN $FCFLAGS"; export FCFLAGS])

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 lithomop3d/Makefile
                 lithomop3d/liblithomop3d/Makefile
                 lithomop3d/module/Makefile])
AC_OUTPUT

dnl end of configure.ac
