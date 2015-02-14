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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = bc
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = _bcmodule.la

subpkgpyexec_PYTHON = bc.py

swig_sources = \
	bc.i \
	BoundaryCondition.i \
	BoundaryConditionPoints.i \
	BCIntegratorSubMesh.i \
	TimeDependent.i \
	TimeDependentPoints.i \
	DirichletBC.i \
	DirichletBoundary.i \
	AbsorbingDampers.i \
	Neumann.i \
	PointForce.i \
	../feassemble/Constraint.i \
	../feassemble/Quadrature.i \
	../feassemble/Integrator.i

swig_generated = \
	bc_wrap.cxx \
	bc.py

_bcmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__bcmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_bcmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _bcmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_bcmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/bc_wrap.cxx $(srcdir)/bc.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/bc_wrap.cxx \
	$(srcdir)/bc.py


# End of file 
