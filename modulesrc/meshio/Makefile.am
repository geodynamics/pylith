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

subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _meshiomodule.la

subpkgpyexec_PYTHON = meshio.py

swig_sources = \
	meshio.i \
	MeshIOObj.i \
	MeshIOAscii.i \
	MeshIOLagrit.i

swig_generated = \
	meshio_wrap.cxx \
	meshio.py

_meshiomodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__meshiomodule_la_SOURCES = $(swig_sources) $(swig_generated)

_meshiomodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _meshiomodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_meshiomodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/meshio_wrap.cxx $(srcdir)/meshio.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/meshio_wrap.cxx \
	$(srcdir)/meshio.py


# End of file 
