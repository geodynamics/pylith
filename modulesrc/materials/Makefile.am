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

subpackage = materials
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _materialsmodule.la

subpkgpyexec_PYTHON = materials.py

swig_sources = \
	materials.i \
	Material.i \
	ElasticMaterial.i \
	ElasticStrain1D.i

swig_generated = \
	materials_wrap.cxx \
	materials.py

_materialsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__materialsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_materialsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _materialsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_materialsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/materials_wrap.cxx $(srcdir)/materials.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/materials_wrap.cxx \
	$(srcdir)/materials.py


# End of file 
