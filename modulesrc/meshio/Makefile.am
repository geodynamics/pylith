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

subpackage = meshio
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpyexec_LTLIBRARIES = _meshiomodule.la

subpkgpyexec_PYTHON = meshio.py

swig_sources = \
	meshio.i \
	MeshIOObj.i \
	MeshIOAscii.i \
	MeshIOLagrit.i \
	MeshIOCubit.i \
	VertexFilter.i \
	VertexFilterVecNorm.i \
	CellFilter.i \
	CellFilterAvg.i \
	DataWriter.i \
	DataWriterVTK.i \
	OutputManager.i \
	OutputSolnSubset.i \
	OutputSolnPoints.i

swig_generated = \
	meshio_wrap.cxx \
	meshio.py

_meshiomodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__meshiomodule_la_SOURCES = $(swig_sources) $(swig_generated)

_meshiomodule_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB)
if ENABLE_CUBIT
  _meshiomodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_meshiomodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

if ENABLE_SWIG
$(srcdir)/meshio_wrap.cxx $(srcdir)/meshio.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) $(PYLITH_SWIG_CPPFLAGS) -Wall -c++ -python $<
endif

if ENABLE_HDF5
  swig_sources += \
	DataWriterHDF5.i \
	DataWriterHDF5Ext.i \
	Xdmf.i
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/meshio_wrap.cxx \
	$(srcdir)/meshio.py


# End of file 
