# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = faults
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _faults.la

subpkgpython_PYTHON = faults.py

swig_sources = \
	faults.i \
	../utils/PyreComponent.i \
	../problems/Physics.i \
	FaultCohesive.i \
	FaultCohesiveKin.i \
	FaultCohesiveImpulses.i \
	KinSrc.i \
	KinSrcStep.i \
	KinSrcRamp.i \
	KinSrcConstRate.i \
	KinSrcBrune.i \
	KinSrcLiuCos.i \
	KinSrcTimeHistory.i

swig_generated = \
	faults_wrap.cxx \
	faults.py

_faults_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__faults_la_SOURCES = $(swig_sources) $(swig_generated)

_faults_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _faults_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/faults_wrap.cxx $(srcdir)/faults.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
else
$(srcdir)/faults_wrap.cxx $(srcdir)/faults.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/faults_wrap.cxx \
	$(srcdir)/faults.py


# End of file
