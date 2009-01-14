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
#

AC_PREREQ(2.59)
AC_INIT([swigtest], [0.1.0], [baagaard@usgs.gov])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_AUX_DIR([./aux-config])
AM_INIT_AUTOMAKE([foreign])

# ----------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_CC

AC_PROG_LIBTOOL
if test "$allow_undefined_flag" = unsupported; then
    # See issue119.
    AM_LDFLAGS="-no-undefined $AM_LDFLAGS"
fi
AM_CONDITIONAL([NO_UNDEFINED], [test "$allow_undefined_flag" = unsupported])
AC_SUBST(AM_LDFLAGS)

# ----------------------------------------------------------------------
# PYTHIA
AC_LANG(C++)

# Check for Python modules and packages.
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])
CIT_PYTHON_SYSCONFIG
AM_PATH_PYTHON([2.3])

# NUMPY
CIT_NUMPY_PYTHON_MODULE
CIT_NUMPY_INCDIR


# SWIG
AC_PROG_SWIG

# ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
	libsrc/Makefile
	modulesrc/Makefile
	python/Makefile
	applications/Makefile])

AC_OUTPUT


# End of file 
