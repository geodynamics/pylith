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

pkgpyexec_LTLIBRARIES = _geometrymodule.la

_geometrymodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

pkgpyexec_PYTHON = geometry.py

swig_sources = \
	geometry.i \
	Bar.i \
	Shape.i \
	Sphere.i \
	Scene.i

swig_generated = \
	geometry.py \
	geometry_wrap.h \
	geometry_wrap.cxx

dist__geometrymodule_la_SOURCES = $(swig_sources) $(swig_generated)

_geometrymodule_la_LIBADD = \
	$(top_builddir)/libsrc/libswigtest.la
if NO_UNDEFINED
_geometrymodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES = -I$(top_srcdir)/libsrc -I$(PYTHON_INCDIR)

$(srcdir)/geometry_wrap.cxx $(srcdir)/geometry.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/geometry_wrap.h \
	$(srcdir)/geometry_wrap.cxx \
	$(srcdir)/geometry.py


# End of file 
