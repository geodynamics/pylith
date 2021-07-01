# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	DataWriter.hh \
	HDF5.hh \
	Xdmf.hh \
	DataWriterHDF5.hh \
	DataWriterHDF5.icc \
	DataWriterHDF5Ext.hh \
	DataWriterHDF5Ext.icc \
	DataWriterVTK.hh \
	DataWriterVTK.icc \
	MeshBuilder.hh \
	MeshIO.hh \
	MeshIO.icc \
	MeshIOAscii.hh \
	MeshIOAscii.icc \
	MeshIOLagrit.hh \
	MeshIOLagrit.icc \
	OutputObserver.hh \
	OutputSubfield.hh \
	OutputSoln.hh \
	OutputSolnDomain.hh \
	OutputSolnBoundary.hh \
	OutputSolnPoints.hh \
	OutputPhysics.hh \
	OutputTrigger.hh \
	OutputTriggerStep.hh \
	OutputTriggerTime.hh \
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
	PsetFileBinary.hh \
	ExodusII.hh


# End of file
