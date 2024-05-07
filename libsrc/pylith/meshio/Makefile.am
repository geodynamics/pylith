# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


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
	MeshConverter.hh \
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
