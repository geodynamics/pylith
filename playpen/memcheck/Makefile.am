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

ACLOCAL_AMFLAGS = -I ./m4

bin_PROGRAMS = runtet4

# Primary source files
runtet4_SOURCES = \
	runtet4.cc

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

runtet4_LDADD = \
	-lpylith \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

CLEANFILES = $(bin_PROGRAMS)

# End of file 
