# -*- makefile -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


subpackage = meshio
include $(top_srcdir)/subpackage.am

subpkginclude_HEADERS = \
	DataWriter.hh \
	HDF5.hh \
	HDF5.icc \
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
	MeshIOCubit.hh
endif

dist_noinst_HEADERS = \
	BinaryIO.hh \
	ExodusII.hh


# End of file
