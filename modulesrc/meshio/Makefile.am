# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _meshio.la

subpkgpython_PYTHON = meshio.py

swig_sources = \
	meshio.i \
	MeshIOObj.i \
	MeshIOAscii.i \
	MeshIOLagrit.i \
	MeshIOCubit.i \
	OutputTrigger.i \
	OutputTriggerStep.i \
	OutputTriggerTime.i \
	DataWriter.i \
	DataWriterHDF5.i \
	DataWriterHDF5Ext.i \
	DataWriterVTK.i \
	OutputObserver.i \
	OutputSoln.i \
	OutputSolnDomain.i \
	OutputSolnBoundary.i \
	OutputSolnPoints.i \
	../utils/PyreComponent.i \
	../problems/ObserverSoln.i \
	OutputPhysics.i \
	../problems/ObserverPhysics.i

swig_generated = \
	meshio_wrap.cxx \
	meshio.py

_meshio_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__meshio_la_SOURCES = $(swig_sources) $(swig_generated)

_meshio_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
if ENABLE_CUBIT
  _meshio_la_LIBADD += -lnetcdf
endif

if ENABLE_SWIG
$(srcdir)/meshio_wrap.cxx $(srcdir)/meshio.py: $(swig_sources)
	$(SWIG) $(PETSC_CC_INCLUDES) $(PYLITH_SWIG_CPPFLAGS) -Wall -c++ -python $<
else
$(srcdir)/meshio_wrap.cxx $(srcdir)/meshio.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/meshio_wrap.cxx \
	$(srcdir)/meshio.py


# End of file
