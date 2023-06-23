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

include $(top_srcdir)/tests/check_cppunit.am

SUBDIRS = data

TESTS = test_meshio test_vtk test_hdf5

check_PROGRAMS = $(TESTS)

# general meshio
test_meshio_SOURCES = \
	TestMeshIO.cc \
	TestMeshIOAscii.cc \
	TestMeshIOAscii_Cases.cc \
	TestMeshIOPetsc.cc \
	TestMeshIOPetsc_Cases.cc \
	TestMeshIOLagrit.cc \
	TestMeshIOLagrit_Cases.cc \
	TestOutputTriggerStep.cc \
	TestOutputTriggerTime.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc
	

# VTK data writer
test_vtk_SOURCES = \
	FieldFactory.cc \
	TestDataWriter.cc \
	TestDataWriterMesh.cc \
	TestDataWriterMaterial.cc \
	TestDataWriterSubmesh.cc \
	TestDataWriterPoints.cc \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMesh_Cases.cc \
	TestDataWriterVTKMaterial.cc \
	TestDataWriterVTKMaterial_Cases.cc \
	TestDataWriterVTKSubmesh.cc \
	TestDataWriterVTKSubmesh_Cases.cc \
	TestDataWriterVTKPoints.cc \
	TestDataWriterVTKPoints_Cases.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc


# :TODO: @brad
# TestDataWriterFaultMesh.cc \
# TestDataWriterVTKFaultMesh.cc \
# TestDataWriterVTKFaultMesh_Cases.cc \
# TestOutputObserver.cc \
# TestOutputSolnDomain.cc \
# TestOutputSolnBoundary.cc \
# TestOutputSolnPoints.cc


test_hdf5_SOURCES = \
	FieldFactory.cc \
	TestDataWriter.cc \
	TestDataWriterMesh.cc \
	TestDataWriterMaterial.cc \
	TestDataWriterSubmesh.cc \
	TestDataWriterPoints.cc \
	TestHDF5.cc \
	TestDataWriterHDF5.cc \
	TestDataWriterHDF5Mesh.cc \
	TestDataWriterHDF5Mesh_Cases.cc \
	TestDataWriterHDF5Material.cc \
	TestDataWriterHDF5Material_Cases.cc \
	TestDataWriterHDF5Submesh.cc \
	TestDataWriterHDF5Submesh_Cases.cc \
	TestDataWriterHDF5Points.cc \
	TestDataWriterHDF5Points_Cases.cc \
	TestDataWriterHDF5ExtMesh.cc \
	TestDataWriterHDF5ExtMesh_Cases.cc \
	TestDataWriterHDF5ExtMaterial.cc \
	TestDataWriterHDF5ExtMaterial_Cases.cc \
	TestDataWriterHDF5ExtSubmesh.cc \
	TestDataWriterHDF5ExtSubmesh_Cases.cc \
	TestDataWriterHDF5ExtPoints.cc \
	TestDataWriterHDF5ExtPoints_Cases.cc \
	$(top_srcdir)/tests/src/FaultCohesiveStub.cc \
	$(top_srcdir)/tests/src/StubMethodTracker.cc \
	$(top_srcdir)/tests/src/driver_cppunit.cc


# TestDataWriterHDF5FaultMesh.cc \
# TestDataWriterHDF5FaultMeshCases.cc \
# TestDataWriterHDF5ExtFaultMesh.cc \
# TestDataWriterHDF5ExtFaultMeshCases.cc


dist_noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOPetsc.hh \
	TestOutputTriggerStep.hh \
	TestOutputTriggerTime.hh \
	FieldFactory.hh \
	TestOutputManager.hh \
	TestOutputSolnSubset.hh \
	TestOutputSolnPoints.hh \
	TestVertexFilterVecNorm.hh \
	TestDataWriter.hh \
	TestDataWriterMesh.hh \
	TestDataWriterMaterial.hh \
	TestDataWriterSubmesh.hh \
	TestDataWriterPoints.hh \
	TestDataWriterVTK.hh \
	TestDataWriterVTKMesh.hh \
	TestDataWriterVTKMaterial.hh \
	TestDataWriterVTKSubmesh.hh \
	TestDataWriterVTKPoints.hh \
	TestHDF5.hh \
	TestDataWriterHDF5.hh \
	TestDataWriterHDF5Mesh.hh \
	TestDataWriterHDF5Material.hh \
	TestDataWriterHDF5Submesh.hh \
	TestDataWriterHDF5Points.hh \
	TestDataWriterHDF5ExtMesh.hh \
	TestDataWriterHDF5ExtMaterial.hh \
	TestDataWriterHDF5ExtSubmesh.hh \
	TestDataWriterHDF5ExtPoints.hh

# :TODO: @brad
# TestDataWriterFaultMesh.hh \
# TestDataWriterVTKFaultMesh.hh
# TestDataWriterHDF5FaultMesh.hh \
# TestDataWriterHDF5ExtFaultMesh.hh


if ENABLE_CUBIT
  test_meshio_SOURCES += \
	TestExodusII.cc \
	TestMeshIOCubit.cc \
	TestMeshIOCubit_Cases.cc
  dist_noinst_HEADERS += \
	TestExodusII.hh \
	TestMeshIOCubit.hh
endif


noinst_TMP = 

#include $(top_srcdir)/tests/data.am

clean-local:
	$(RM) $(RM_FLAGS) mesh*.txt *.h5 *.xmf *.dat *.dat.info *.vtk


# End of file
