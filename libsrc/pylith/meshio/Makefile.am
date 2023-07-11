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
# Copyright (c) 2010-2022 University of California, Davis
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
	MeshIOAscii.hh \
	MeshIOAscii.icc \
	MeshIOPetsc.hh \
	MeshIOPetsc.icc \
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

dist_noinst_HEADERS = \
	BinaryIO.hh \
	ExodusII.hh


# End of file
