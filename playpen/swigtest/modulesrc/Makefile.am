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

swig_sources = \
	geometry.i \
	Bar.i \
	Shape.i \
	Sphere.i \
	Scene.i

dist__geometrymodule_la_SOURCES = $(swig_sources)

pkgpyexec_PYTHON = geometry.py

_geometrymodule_la_SOURCES = \
	geometry_wrap.cxx

_geometrymodule_la_LIBADD = \
	$(top_builddir)/libsrc/libswigtest.la
if NO_UNDEFINED
_geometrymodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES = -I$(top_srcdir)/libsrc -I$(PYTHON_INCDIR)

$(srcdir)/geometry_wrap.cxx $(srcdir)/geometry.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


CLEANFILES = $(srcdir)/geometry_wrap.cxx $(srcdir)/geometry.py


# End of file 
