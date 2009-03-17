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

subpkginclude_HEADERS = \
	Formulation.hh \
	Solver.hh \
	SolverLinear.hh \
	SolverNonlinear.hh \
	problemsfwd.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
