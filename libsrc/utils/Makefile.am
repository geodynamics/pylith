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
	Endian.hh \
	array.hh \
	arrayfwd.hh \
	petscfwd.h \
	sievefwd.hh \
	sievetypes.hh

noinst_HEADERS = 

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS

# End of file 
