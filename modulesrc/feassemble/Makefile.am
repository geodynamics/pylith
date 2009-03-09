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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _feassemblemodule.la

subpkgpyexec_PYTHON = feassemble.py

swig_sources = \
	feassemble.i \
	CellGeometry.i \
	GeometryPoint1D.i

swig_generated = \
	feassemble_wrap.cxx \
	feassemble.py

_feassemblemodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__feassemblemodule_la_SOURCES = $(swig_sources) $(swig_generated)

_feassemblemodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _feassemblemodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_feassemblemodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/feassemble_wrap.cxx $(srcdir)/feassemble.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/feassemble_wrap.cxx \
	$(srcdir)/feassemble.py


# End of file 
