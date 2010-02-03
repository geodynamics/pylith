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

subpackage = friction
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _frictionmodule.la

subpkgpyexec_PYTHON = friction.py

swig_sources = \
	friction.i \
	FrictionModel.i \
	StaticFriction.i

swigincludedir = $(pkgdatadir)/swig/$(subpackage)
swiginclude_HEADERS = \
	FrictionModel.i

swig_generated = \
	friction_wrap.cxx \
	friction.py

_frictionmodule_la_LDFLAGS = -module -avoid-version \
	-g $(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__frictionmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_frictionmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _frictionmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_frictionmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += -g $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

if ENABLE_SWIG
$(srcdir)/friction_wrap.cxx $(srcdir)/friction.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/friction_wrap.cxx \
	$(srcdir)/friction.py


# End of file 
