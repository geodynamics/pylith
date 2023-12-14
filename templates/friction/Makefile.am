# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
#

ACLOCAL_AMFLAGS = -I ./m4

SUBDIRS = tests

subpkgpythondir = $(pythondir)/pylith/friction/contrib


# LIBRARY --------------------------------------------------------------
lib_LTLIBRARIES = libfrictioncontrib.la

libfrictioncontrib_la_SOURCES = \
	ViscousFriction.cc

dist_noinst_HEADERS = \
	ViscousFriction.hh

libfrictioncontrib_la_LDFLAGS = $(AM_LDFLAGS)

libfrictioncontrib_la_LIBADD = \
	-lpylith \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS = $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) 
AM_CPPFLAGS += $(PETSC_CC_INCLUDES)

# MODULE ---------------------------------------------------------------

subpkgpython_LTLIBRARIES = _frictioncontribmodule.la

subpkgpython_PYTHON = \
	frictioncontrib.py \
	__init__.py

swig_sources = \
	frictioncontrib.i \
	ViscousFriction.i

swig_generated = \
	frictioncontrib_wrap.cxx \
	frictioncontrib.py

_frictioncontribmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__frictioncontribmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_frictioncontribmodule_la_LIBADD = libfrictioncontrib.la \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

AM_CPPFLAGS += -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR)

$(srcdir)/frictioncontrib_wrap.cxx $(srcdir)/frictioncontrib.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $(SWIG_FLAGS) $<


MAINTAINERCLEANFILES = $(swig_generated)



# PYTHON ---------------------------------------------------------------

nobase_subpkgpython_PYTHON = \
	__init__.py \
	ViscousFriction.py


# End of file 
