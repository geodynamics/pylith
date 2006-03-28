#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $

AC_PREREQ(2.59)
AC_INIT([Lithomop3d], [0.7.1], [willic3@rpi.edu])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_SRCDIR([lithomop3d/applications/lithomop3dapp.py])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# 'configure' options
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
AC_ARG_WITH([embedding],
    [AC_HELP_STRING([--with-embedding],
        [embed Python with Lithomop in a single executable @<:@default=yes@:>@])],
    [want_embedding="$withval"],
    [want_embedding=yes])
AM_CONDITIONAL([COND_EMBEDDING], [test "$want_embedding" = yes])

# Checks for programs.
AM_PATH_PYTHON([2.3])
CIT_PYTHON_SYSCONFIG
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
PYTHONPATH="${pyexecdir}${PYTHONPATH:+:${PYTHONPATH}}"; export PYTHONPATH
AC_SUBST([PYTHONPATH])

# Set the language for MPI tests.
AC_LANG(C++)

# Checks for libraries.
CIT_CHECK_LIB_MPI

# Checks for header files.
CIT_HEADER_MPI

# Checks for typedefs, structures, and compiler characteristics.
AC_LANG(Fortran)
fc_append=no
AC_MSG_CHECKING([whether $FC supports OPEN control item 'position="append"'])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="old",position="append")]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_POSITION_APPEND $FCFLAGS"; export FCFLAGS
    fc_append=yes
], [
    AC_MSG_RESULT(no)
])
AC_MSG_CHECKING([whether $FC supports OPEN control item 'access="append"'])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="old",access="append")]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_ACCESS_APPEND $FCFLAGS"; export FCFLAGS
    fc_append=yes
], [
    AC_MSG_RESULT(no)
])
AS_IF([test $fc_append = yes], [], [
    AC_MSG_FAILURE([cannot determine method for appending to Fortran files])
])
AC_MSG_CHECKING([whether $FC supports stream i/o])
AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="stream",
     & form="unformatted")
      write(10,pos=1) 1,2,3.0d0]])
], [
    AC_MSG_RESULT(yes)
    FCFLAGS="-DFORTRAN_STREAM_IO $FCFLAGS"; export FCFLAGS
], [
        AC_MSG_RESULT(no)
        AC_MSG_CHECKING([whether $FC supports f77-style binary direct-access i/o])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="direct",recl=1,
     & form="unformatted")
      write(10,rec=1) 1,2,3.0d0]])
    ], [
        AC_MSG_RESULT(yes)
        FCFLAGS="-DFORTRAN_F77_IO $FCFLAGS"; export FCFLAGS
        AC_MSG_CHECKING([whether $FC supports I/O specifiers 'advance' and 'eor'])
        AC_COMPILE_IFELSE([
            AC_LANG_PROGRAM([], [[      open(10,file="foo",status="new",access="direct",recl=1,
     & form="unformatted")
      write(10,rec=1,advance='yes',eor=10) 1,2,3.0d0
 10   continue]])
        ], [
            AC_MSG_RESULT(yes)
            FCFLAGS="-DFORTRAN_EOR $FCFLAGS"; export FCFLAGS
        ], [
            AC_MSG_RESULT(no)
        ])
    ], [
        AC_MSG_RESULT(no)
        AC_MSG_WARN([cannot determine how to produce binary direct-access files with variable record length])
        FCFLAGS="-DFORTRAN_NO_BINARY $FCFLAGS"; export FCFLAGS
    ])
])

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 lithomop3d/Makefile
                 lithomop3d/liblithomop3d/Makefile
                 lithomop3d/module/Makefile])
AC_OUTPUT

dnl end of configure.ac
