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

subpackage = mpi
include $(top_srcdir)/subpackage.am

subpkgpyexec_LTLIBRARIES = _mpimodule.la

subpkgpyexec_PYTHON = mpi.py

swig_sources = \
	mpi.i \
	mpi_comm.i

swig_generated = \
	mpi_wrap.cxx \
	mpi.py

_mpimodule_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__mpimodule_la_SOURCES = $(swig_sources) $(swig_generated)

_mpimodule_la_LIBADD = \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIBS)
if ENABLE_CUBIT
  _mpimodule_la_LIBADD += -lnetcdf_c++ -lnetcdf
endif
if NO_UNDEFINED
_mpimodule_la_LIBADD += \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)
endif

INCLUDES += $(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) $(PETSC_INCLUDE)

$(srcdir)/mpi_wrap.cxx $(srcdir)/mpi.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<


MAINTAINERCLEANFILES = \
	$(srcdir)/mpi_wrap.cxx \
	$(srcdir)/mpi.py


# End of file 
