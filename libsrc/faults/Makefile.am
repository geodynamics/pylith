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

subpackage = faults
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	BruneSlipFn.hh \
	BruneSlipFn.icc \
	CohesiveTopology.hh \
	EqKinSrc.hh \
	Fault.hh \
	Fault.icc \
	FaultCohesive.hh \
	FaultCohesiveKin.hh \
	FaultCohesiveKin.icc \
	SlipTimeFn.hh

noinst_HEADERS =

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
