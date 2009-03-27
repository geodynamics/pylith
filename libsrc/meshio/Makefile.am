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

subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	CellFilter.hh \
	CellFilterAvg.hh \
	DataWriter.hh \
	DataWriter.cc \
	DataWriterVTK.hh \
	DataWriterVTK.icc \
	DataWriterVTK.cc \
	MeshBuilder.hh \
	MeshIO.hh \
	MeshIO.icc \
	MeshIOAscii.hh \
	MeshIOAscii.icc \
	MeshIOLagrit.hh \
	MeshIOLagrit.icc \
	OutputManager.hh \
	OutputManager.cc \
	OutputSolnSubset.hh \
	UCDFaultFile.hh \
	VertexFilter.hh \
	VertexFilterVecNorm.hh \
	meshiofwd.hh

if ENABLE_CUBIT
  subpkginclude_HEADERS += \
	MeshIOCubit.hh \
	MeshIOCubit.icc
endif

noinst_HEADERS = \
	BinaryIO.hh \
	GMVFile.hh \
	GMVFileAscii.hh \
	GMVFileAscii.icc \
	GMVFileBinary.hh \
	PsetFile.hh \
	PsetFileAscii.hh \
	PsetFileAscii.icc \
	PsetFileBinary.hh

# export
clean-local: clean-subpkgincludeHEADERS
BUILT_SOURCES = export-subpkgincludeHEADERS
CLEANFILES = export-subpkgincludeHEADERS


# End of file 
