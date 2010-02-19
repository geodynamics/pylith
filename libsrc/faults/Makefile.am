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
	FaultCohesive.icc \
	FaultCohesiveLagrange.hh \
	FaultCohesiveTract.hh \
	FaultCohesiveDyn.hh \
	FaultCohesiveKin.hh \
	LiuCosSlipFn.hh \
	LiuCosSlipFn.icc \
	SlipTimeFn.hh \
	StepSlipFn.hh \
	StepSlipFn.icc \
	TimeHistorySlipFn.hh \
	TimeHistorySlipFn.icc \
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
