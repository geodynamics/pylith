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

subpackage = problems
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _problemsmodule.la

subpkgpyexec_PYTHON = problems.py

swig_sources = \
	problems.i \
	Formulation.i \
	Solver.i \
	SolverLinear.i \
	SolverNonlinear.i


swig_generated = \
	problems_wrap.cxx \
	problems.py

_problemsmodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__problemsmodule_la_SOURCES = $(swig_sources) $(swig_generated)

_problemsmodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _problemsmodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_problemsmodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(NUMPY_INCDIR) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/problems_wrap.cxx $(srcdir)/problems.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/problems_wrap.cxx \
	$(srcdir)/problems.py


# End of file 
