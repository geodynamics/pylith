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

ACLOCAL_AMFLAGS = -I ./m4

SUBDIRS = tests

subpkgpythondir = $(pythondir)/pylith/materials/contrib


# LIBRARY --------------------------------------------------------------
lib_LTLIBRARIES = libmaterialscontrib.la

libmaterialscontrib_la_SOURCES = \
	PlaneStrainState.cc

noinst_HEADERS = \
	PlaneStrainState.hh

libmaterialscontrib_la_LDFLAGS = $(AM_LDFLAGS)

libmaterialscontrib_la_LIBADD = \
	-lpylith \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) 
AM_CPPFLAGS += $(PETSC_CC_INCLUDES)

# MODULE ---------------------------------------------------------------

subpkgpython_LTLIBRARIES = _materialscontribmodule.la

subpkgpython_PYTHON = \
	materialscontrib.py \
	__init__.py

swig_sources = \
	materialscontrib.i \
	PlaneStrainState.i

swig_generated = \
	materialscontrib_wrap.cxx \
	materialscontrib.py

_materialscontribmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__materialscontribmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_materialscontribmodule_la_LIBADD = libmaterialscontrib.la \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS += -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR)

$(srcdir)/materialscontrib_wrap.cxx $(srcdir)/materialscontrib.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $(SWIG_FLAGS) $<


MAINTAINERCLEANFILES = $(swig_generated)



# PYTHON ---------------------------------------------------------------

nobase_subpkgpython_PYTHON = \
	__init__.py \
	PlaneStrainState.py


# End of file 
