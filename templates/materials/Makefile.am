# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
#

ACLOCAL_AMFLAGS = -I ./m4

SUBDIRS = tests

subpkgpythondir = $(pythondir)/pylith/materials/contrib


# LIBRARY --------------------------------------------------------------
lib_LTLIBRARIES = libmaterialscontrib.la

libmaterialscontrib_la_SOURCES = \
	PlaneStrainState.cc

dist_noinst_HEADERS = \
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
