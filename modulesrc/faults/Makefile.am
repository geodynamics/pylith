# -*- Makefile -*-
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
	EqKinSrc.i \
	Fault.i \
	FaultCohesive.i \
	FaultCohesiveDyn.i \
	FaultCohesiveKin.i \
	../topology/SubMesh.i \
	../feassemble/Quadrature.i \
	../feassemble/Integrator.i

swig_generated = \
	faults_wrap.cxx \
	faults.py

_faultsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__faultsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_faultsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _faultsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_faultsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/faults_wrap.cxx $(srcdir)/faults.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/faults_wrap.cxx \
	$(srcdir)/faults.py


# End of file 
