#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.1 2005/08/24 22:45:51 leif Exp $

AC_PREREQ(2.59)
AC_INIT([PyLith3D], [0.8.2], [cig-short@geodynamics.org])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_SRCDIR([pylith3d/applications/pylith3dapp.py])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])

# 'configure' options
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
AC_ARG_WITH([embedding],
    [AC_HELP_STRING([--with-embedding],
        [embed Python with PyLith in a single executable @<:@default=yes@:>@])],
    [want_embedding="$withval"],
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

# Checks for Python modules and packages.

builddir=`pwd`
save_PYTHONPATH="$PYTHONPATH"
PYTHONPATH="$builddir/python:$PYTHONPATH"; export PYTHONPATH
cd $srcdir

AC_MSG_NOTICE([downloading missing Python dependencies])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py install_deps -zmxd $builddir/deps >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_FAILURE([cannot download missing Python dependencies])])

AC_MSG_NOTICE([building Python dependencies])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py develop -H None -f $builddir/deps -x -d $builddir/python >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [],
      [AC_MSG_FAILURE([building Python dependencies])])

AC_MSG_CHECKING([for egg-related flags])
AS_IF([AC_TRY_COMMAND([$PYTHON setup.py egg_flags >&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD])],
      [AC_MSG_RESULT(ok)
       . egg-flags.sh
       rm -f egg-flags.sh
      ],
      [AC_MSG_RESULT(failed)
      AC_MSG_FAILURE([cannot scan Python eggs for flags])])

cd $builddir
PYTHONPATH="$save_PYTHONPATH"
PYTHONPATH="${pythondir}:${pyexecdir}${save_PYTHONPATH:+:${save_PYTHONPATH}}"

AC_SUBST(PYTHONPATH)
AC_SUBST(PYTHON_EGG_CFLAGS)
AC_SUBST(PYTHON_EGG_CPPFLAGS)
AC_SUBST(PYTHON_EGG_LDFLAGS)

# Checks for programs.
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
F77=$FC
FFLAGS=$FCFLAGS
AC_PROVIDE([AC_PROG_F77])
AC_FC_WRAPPERS
AC_PROG_LIBTOOL

# Set the language for MPI tests.
AC_LANG(C++)

# Checks for libraries.
CIT_CHECK_LIB_MPI
CIT_CHECK_LIB_PETSC
CIT_CHECK_LIB_PETSC_SIEVE

# Checks for header files.
CIT_HEADER_MPI
CIT_HEADER_PETSC

# Checks for typedefs, structures, and compiler characteristics.
CIT_FC_OPEN_APPEND
CIT_FC_STREAM_IO

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 pylith3d/Makefile
                 pylith3d/libpylith3d/Makefile
                 pylith3d/module/Makefile
		 pylith3d/utils/Makefile])
AC_OUTPUT

dnl end of configure.ac
