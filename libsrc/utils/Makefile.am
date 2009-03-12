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

subpackage = utils
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	EventLogger.hh \
	EventLogger.icc \
	TestArray.hh \
	array.hh \
	arrayfwd.hh \
	constdefs.h \
	lapack.h \
	macrodefs.h \
	petscfwd.h \
	sievefwd.hh \
	sievetypes.hh \
	vectorfields.hh \
	utilsfwd.hh

noinst_HEADERS = 

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# End of file 
