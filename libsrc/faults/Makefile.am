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
	ConstRateSlipFn.hh \
	ConstRateSlipFn.icc \
	CohesiveTopology.hh \
	EqKinSrc.hh \
	Fault.hh \
	Fault.icc \
	FaultCohesive.hh \
	FaultCohesiveDyn.hh \
	FaultCohesiveDyn.icc \
	FaultCohesiveKin.hh \
	FaultCohesiveKin.icc \
	LiuCosSlipFn.hh \
	LiuCosSlipFn.icc \
	SlipTimeFn.hh \
	StepSlipFn.hh \
	StepSlipFn.icc \
	faultsfwd.hh

noinst_HEADERS = \
	TopologyOps.hh \
	TopologyVisitors.hh \
	TopologyVisitors.cc

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
