dnl ================================================================================================
dnl This code is part of PyLith, developed through the Computational Infrastructure
dnl for Geodynamics (https://github.com/geodynamics/pylith).
dnl
dnl Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
dnl All rights reserved.
dnl
dnl See https://mit-license.org/ and LICENSE.md and for license information. 
dnl ================================================================================================
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([PyLith], [5.0.0dev], [https://geodynamics.org/resources/pylith])
AC_CONFIG_AUX_DIR([./aux-config])
AC_CONFIG_HEADER([portinfo])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects tar-pax])

dnl ----------------------------------------------------------------------
dnl PYTHON
AC_ARG_VAR(PYTHON, [Python interpreter])
AC_ARG_VAR(PYTHONPATH, [Python module search path])

dnl SWIG
AC_ARG_ENABLE([swig],
    [AC_HELP_STRING([--enable-swig],
        [enable generating modules with SWIG) @<:@default=no@:>@])],
	[if test "$enableval" = yes ; then enable_swig=yes; else enable_swig=no; fi],
	[enable_swig=no])
AM_CONDITIONAL([ENABLE_SWIG], [test "$enable_swig" = yes])

dnl TESTING w/cppunit and Python unittest
AC_ARG_ENABLE([testing],
	[  --enable-testing        Enable unit testing with CppUnit (requires CppUnit) [[default=no]]],
	[if test "$enableval" = yes ; then enable_testing=yes; else enable_testing=no; fi],
	[enable_testing=no])
AM_CONDITIONAL([ENABLE_TESTING], [test "$enable_testing" = yes])

dnl TEST COVERAGE w/locv and python-coverage
AC_ARG_ENABLE([test-coverage],
	[  --enable-test-coverage  Enable test coverage with lcov and python-coverage [[default=no]]],
	[if test "$enableval" = yes ; then enable_test_coverage=yes; else enable_test_coverage=no; fi],
	[enable_test_coverage=no])
AC_ARG_WITH([python-coverage],
    [AC_HELP_STRING([--with-python-coverage],
        [set executable for python-coverage @<:@default=coverage2@:>@])],
	[python_coverage=$withval],
	[python_coverage="coverage2"])
AC_SUBST(python_coverage)

dnl sanitizer
AC_ARG_ENABLE([sanitizer],
    [AC_HELP_STRING([--enable-sanitizer],
        [enable compiling C++ code with -fsanitizer @<:@default=no@:>@])],
	[if test "$enableval" != no; then enable_sanitizer=yes; else enable_sanitizer=no; fi;
	if test "$enableval" = yes; then sanitizer=address; else sanitizer=$enableval; fi],
	[enable_sanitizer=no])
if test "$enable_sanitizer" = "yes"; then
  SANITIZER_CXXFLAGS="-fsanitize=$sanitizer"; export SANITIZER_CXXFLAGS
fi
AC_SUBST(SANITIZER_CXXFLAGS)

dnl CUBIT I/O w/netcdf
AC_ARG_ENABLE([cubit],
    [AC_HELP_STRING([--enable-cubit],
        [enable reading/writing Cubit EXODUS files (requires netcdf) @<:@default=yes@:>@])],
	[if test "$enableval" = yes; then enable_cubit=yes; else enable_cubit=no; fi],
	[enable_cubit=yes])
if test "$enable_cubit" = "yes"; then
  CPPFLAGS="-DENABLE_CUBIT $CPPFLAGS"; export CPPFLAGS
  PYLITH_SWIG_CPPFLAGS="-DENABLE_CUBIT $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS
fi
AM_CONDITIONAL([ENABLE_CUBIT], [test "$enable_cubit" = yes])
AC_SUBST(PYLITH_SWIG_CPPFLAGS)

dnl PETSc w/HDF5
AC_ARG_ENABLE([hdf5],
    [AC_HELP_STRING([--enable-hdf5],
        [enable output to HDF5 files via PETSc @<:@default=yes@:>@])],
	[if test "$enableval" = yes; then enable_hdf5=yes; else enable_hdf5=no; fi],
	[enable_hdf5=yes])
if test "$enable_hdf5" = "yes"; then
  CPPFLAGS="-DENABLE_HDF5 $CPPFLAGS"; export CPPFLAGS
  PYLITH_SWIG_CPPFLAGS="-DENABLE_HDF5 $PYLITH_SWIG_CPPFLAGS"; export PYLITH_SWIG_CPPFLAGS
fi
AM_CONDITIONAL([ENABLE_HDF5], [test "$enable_hdf5" = yes])
AC_SUBST(PYLITH_SWIG_CPPFLAGS)


dnl ----------------------------------------------------------------------
dnl C/C++/libtool/install
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CXXCPP
AC_DISABLE_STATIC

AC_PROG_LIBTOOL
AC_PROG_INSTALL

AX_CXX_COMPILE_STDCXX(14)

dnl PYTHON/PYTHIA (nemesis must be set before python)
CIT_PATH_NEMESIS
AM_PATH_PYTHON([3.8])
CIT_PYTHON_SYSCONFIG
CIT_CHECK_PYTHON_HEADER
CIT_CHECK_PYTHON_SHARED
CIT_PYTHON_MODULE([pythia],[1.1.0])

dnl MPI
AC_LANG(C++)
AC_SEARCH_LIBS([MPI_Init], [mpi mpich], [], [AC_MSG_ERROR([MPI library not found])])
AC_CHECK_HEADER([mpi.h], [], [AC_MSG_ERROR([header 'mpi.h' not found])])

dnl PETSC
AC_LANG(C)
CIT_PATH_PETSC([3.23.4])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

dnl Large file support
AC_SYS_LARGEFILE

dnl NUMPY
CIT_NUMPY_PYTHON_MODULE
CIT_NUMPY_INCDIR
CIT_PYTHON_MODULE([numpy],[1.20.0])

dnl SWIG
if test "$enable_swig" = "yes" ; then
  AC_PROG_SWIG([4.0.0])
fi

dnl CATCH2
if test "$enable_testing" = "yes" ; then
  CIT_CATCH2_HEADER
  CIT_CATCH2_LIB
fi

dnl CUBIT (netcdf)
if test "$enable_cubit" = "yes" ; then
  AX_LIB_NETCDF4()
fi

dnl HDF5
if test "$enable_hdf5" = "yes" ; then
  AC_REQUIRE_CPP
  AC_LANG(C++)
  CIT_HDF5_HEADER
  CIT_HDF5_LIB
  CIT_HDF5_LIB_PARALLEL

  CIT_PYTHON_MODULE([h5py],[3.0.0])
fi

dnl PROJ
AC_REQUIRE_CPP
AC_LANG(C)
CIT_PROJ6_HEADER
CIT_PROJ6_LIB

dnl SPATIALDATA
AC_REQUIRE_CPP
AC_LANG(C++)
CIT_SPATIALDATA_HEADER
CIT_SPATIALDATA_LIB
CIT_PYTHON_MODULE([spatialdata],[3.1.0])

dnl TEST COVERAGE
AM_CONDITIONAL([ENABLE_TEST_COVERAGE], [test "$enable_test_coverage" = yes])
if test "$enable_test_coverage" = "yes" ; then
  AC_PATH_PROG(LCOV, lcov)
  if test -z "$LCOV"; then
    AC_MSG_FAILURE([cannot find 'lcov' program for C++ test coverage.])
  fi
  AC_PATH_PROG(PYTHON_COVERAGE, $python_coverage)
  if test -z "$PYTHON_COVERAGE"; then
    AC_MSG_FAILURE([cannot find '$python_coverage' program for Python test coverage.])
  fi
fi

dnl Check whether compiler sets variable with function names.
CIT_FUNCTIONSTRING

dnl VERSION
CIG_PKG_GIT(PYLITH)
AC_DEFINE_UNQUOTED([PYLITH_VERSION], ["$PACKAGE_VERSION"], [Define PyLith version])
AC_DEFINE_UNQUOTED([PYLITH_DOI], ["10.5281/zenodo.14635926"], [Define PyLith doi])

dnl ENDIANNESS
AC_C_BIGENDIAN

dnl ----------------------------------------------------------------------
AC_CONFIG_FILES([Makefile
		pylith/Makefile
		libsrc/Makefile
		libsrc/pylith/Makefile
		libsrc/pylith/bc/Makefile
		libsrc/pylith/feassemble/Makefile
		libsrc/pylith/fekernels/Makefile
		libsrc/pylith/faults/Makefile
		libsrc/pylith/friction/Makefile
		libsrc/pylith/materials/Makefile
		libsrc/pylith/meshio/Makefile
		libsrc/pylith/problems/Makefile
		libsrc/pylith/topology/Makefile
		libsrc/pylith/utils/Makefile
        libsrc/pylith/testing/Makefile
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
		tests/Makefile
		tests/src/Makefile
		tests/libtests/Makefile
		tests/libtests/bc/Makefile
		tests/libtests/bc/data/Makefile
		tests/libtests/faults/Makefile
		tests/libtests/faults/data/Makefile
		tests/libtests/feassemble/Makefile
		tests/libtests/feassemble/data/Makefile
		tests/libtests/materials/Makefile
		tests/libtests/materials/data/Makefile
		tests/libtests/meshio/Makefile
		tests/libtests/problems/Makefile
		tests/libtests/problems/data/Makefile
		tests/libtests/meshio/data/Makefile
		tests/libtests/topology/Makefile
		tests/libtests/topology/data/Makefile
		tests/libtests/testing/Makefile
		tests/libtests/utils/Makefile
		tests/pytests/Makefile
		tests/mmstests/Makefile
		tests/mmstests/linearelasticity/Makefile
		tests/mmstests/linearelasticity/nofaults-2d/Makefile
		tests/mmstests/linearelasticity/nofaults-3d/Makefile
		tests/mmstests/linearelasticity/faults-2d/Makefile
		tests/mmstests/incompressibleelasticity/Makefile
		tests/mmstests/incompressibleelasticity/nofaults-2d/Makefile
		tests/mmstests/poroelasticity/Makefile
		tests/mmstests/poroelasticity/nofaults-2d/Makefile
		tests/fullscale/Makefile
		tests/fullscale/cornercases/Makefile
		tests/fullscale/cornercases/nofaults-2d/Makefile
		tests/fullscale/cornercases/nofaults-3d/Makefile
		tests/fullscale/cornercases/faults-2d/Makefile
		tests/fullscale/linearelasticity/Makefile
		tests/fullscale/linearelasticity/nofaults-2d/Makefile
		tests/fullscale/linearelasticity/nofaults-3d/Makefile
		tests/fullscale/linearelasticity/faults-2d/Makefile
		tests/fullscale/linearelasticity/faults-3d/Makefile
		tests/fullscale/linearelasticity/faults-3d-buried/Makefile
		tests/fullscale/linearelasticity/greensfns-2d/Makefile
		tests/fullscale/incompressibleelasticity/Makefile
		tests/fullscale/incompressibleelasticity/nofaults-2d/Makefile
		tests/fullscale/incompressibleelasticity/nofaults-3d/Makefile
		tests/fullscale/viscoelasticity/Makefile
		tests/fullscale/viscoelasticity/nofaults-2d/Makefile
		tests/fullscale/viscoelasticity/nofaults-3d/Makefile
		tests/fullscale/poroelasticity/Makefile
		tests/fullscale/poroelasticity/terzaghi/Makefile
		tests/fullscale/poroelasticity/mandel/Makefile
		tests/fullscale/poroelasticity/cryer/Makefile
		tests/fullscale/petsc/Makefile
		tests/fullscale/eqinfo/Makefile
		tests/manual/Makefile
		tests/manual/powerlaw-2d/Makefile
		tests/manual/powerlaw-3d/Makefile
		developer/Makefile
		docs/Makefile
		release-notes/Makefile
		examples/Makefile
		examples/barwaves-2d/Makefile
		examples/box-2d/Makefile
		examples/box-3d/Makefile
		examples/reverse-2d/Makefile
		examples/strikeslip-2d/Makefile
		examples/crustal-strikeslip-2d/Makefile
		examples/crustal-strikeslip-3d/Makefile
		examples/subduction-2d/Makefile
		examples/subduction-3d/Makefile
		examples/magma-2d/Makefile
		examples/poroelastic-outerrise-2d/Makefile
		examples/troubleshooting-2d/Makefile
		examples/meshing-cubit/Makefile
		examples/meshing-cubit/surface-nurbs/Makefile
		examples/meshing-cubit/cellsize/Makefile
		templates/Makefile
		share/Makefile
		])


AC_OUTPUT

dnl end of configure.ac
