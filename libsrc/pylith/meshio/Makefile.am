# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	CellFilter.hh \
	CellFilter.cc \
	CellFilterAvg.hh \
	CellFilterAvg.cc \
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
	MeshIOSieve.hh \
	MeshIOSieve.icc \
	OutputManager.hh \
	OutputManager.cc \
	OutputSolnSubset.hh \
	UCDFaultFile.hh \
	VertexFilter.hh \
	VertexFilter.cc \
	VertexFilterVecNorm.hh \
	VertexFilterVecNorm.cc \
	meshiofwd.hh

if ENABLE_HDF5
  subpkginclude_HEADERS += \
	HDF5.hh \
	Xdmf.hh \
	DataWriterHDF5.hh \
	DataWriterHDF5.icc \
	DataWriterHDF5.cc \
	DataWriterHDF5Ext.hh \
	DataWriterHDF5Ext.icc \
	DataWriterHDF5Ext.cc
endif

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
	PsetFileBinary.hh \
	ExodusII.hh


# End of file 
