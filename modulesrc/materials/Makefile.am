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

subpackage = materials
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = _materialsmodule.la

subpkgpyexec_PYTHON = materials.py

swig_sources = \
	materials.i \
	Material.i \
	ElasticMaterial.i \
	ElasticPlaneStrain.i \
	ElasticPlaneStress.i \
	ElasticIsotropic3D.i \
	MaxwellIsotropic3D.i \
	MaxwellPlaneStrain.i \
	GenMaxwellIsotropic3D.i \
	GenMaxwellPlaneStrain.i \
	GenMaxwellQpQsIsotropic3D.i \
	PowerLaw3D.i \
	PowerLawPlaneStrain.i \
	DruckerPrager3D.i \
	DruckerPragerPlaneStrain.i

swigincludedir = $(pkgdatadir)/swig/$(subpackage)
swiginclude_HEADERS = \
	Material.i \
	ElasticMaterial.i

swig_generated = \
	materials_wrap.cxx \
	materials.py

_materialsmodule_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__materialsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_materialsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _materialsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_materialsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/materials_wrap.cxx $(srcdir)/materials.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/materials_wrap.cxx \
	$(srcdir)/materials.py


# End of file 
