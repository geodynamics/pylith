# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = faults
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _faultsmodule.la

subpkgpyexec_PYTHON = faults.py

swig_sources = \
	faults.i \
	SlipTimeFn.i \
	ConstRateSlipFn.i \
	StepSlipFn.i \
	BruneSlipFn.i \
	LiuCosSlipFn.i \
	TimeHistorySlipFn.i \
	EqKinSrc.i \
	TractPerturbation.i \
	Fault.i \
	FaultCohesive.i \
	FaultCohesiveLagrange.i \
	FaultCohesiveKin.i \
	FaultCohesiveDyn.i \
	FaultCohesiveImpulses.i \
	FaultCohesiveTract.i \
	../topology/Mesh.i \
	../feassemble/Quadrature.i \
	../feassemble/Integrator.i

swig_generated = \
	faults_wrap.cxx \
	faults.py

_faultsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__faultsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_faultsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _faultsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_faultsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

AM_CPPFLAGS += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_CC_INCLUDES)

if ENABLE_SWIG
$(srcdir)/faults_wrap.cxx $(srcdir)/faults.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/faults_wrap.cxx \
	$(srcdir)/faults.py


# End of file 
