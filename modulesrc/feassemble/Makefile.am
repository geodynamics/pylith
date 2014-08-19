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

subpackage = feassemble
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = _feassemblemodule.la

subpkgpyexec_PYTHON = feassemble.py

swig_sources = \
	feassemble.i \
	CellGeometry.i \
	GeometryLine2D.i \
	GeometryLine3D.i \
	GeometryTri2D.i \
	GeometryTri3D.i \
	GeometryQuad2D.i \
	GeometryQuad3D.i \
	GeometryTet3D.i \
	GeometryHex3D.i \
	QuadratureRefCell.i \
	Quadrature.i \
	Integrator.i \
	IntegratorElasticity.i \
	ElasticityImplicit.i \
	ElasticityExplicit.i \
	ElasticityExplicitTri3.i \
	ElasticityExplicitTet4.i \
	IntegratorElasticityLgDeform.i \
	ElasticityImplicitLgDeform.i \
	ElasticityExplicitLgDeform.i \
	ElasticityImplicitCUDA.i

swig_generated = \
	feassemble_wrap.cxx \
	feassemble.py

_feassemblemodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__feassemblemodule_la_SOURCES = $(swig_sources) $(swig_generated)

_feassemblemodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _feassemblemodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_feassemblemodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/feassemble_wrap.cxx $(srcdir)/feassemble.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/feassemble_wrap.cxx \
	$(srcdir)/feassemble.py


# End of file 
