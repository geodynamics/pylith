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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmeshio

check_PROGRAMS = testmeshio

# Primary source files
testmeshio_SOURCES = \
	TestMeshIO.cc \
	TestMeshIOAscii.cc \
	TestMeshIOAscii_Cases.cc \
	TestMeshIOLagrit.cc \
	TestMeshIOLagrit_Cases.cc \
	TestOutputTriggerStep.cc \
	TestOutputTriggerTime.cc \
	FieldFactory.cc \
	TestDataWriter.cc \
	TestDataWriterMesh.cc \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMesh_Cases.cc \
	TestDataWriterMaterial.cc \
	TestDataWriterVTKMaterial.cc \
	TestDataWriterVTKMaterial_Cases.cc \
	TestDataWriterSubmesh.cc \
	TestDataWriterVTKSubmesh.cc \
	TestDataWriterVTKSubmesh_Cases.cc \
	test_meshio.cc


#TestCellFilterAvg.cc \
#TestVertexFilterVecNorm.cc

#TestDataWriterSubmesh.cc \
#TestDataWriterBCMesh.cc \
# TestDataWriterVTKSubmesh.cc \
# TestDataWriterVTKSubmeshCases.cc \
# TestDataWriterVTKBCMesh.cc \
# TestDataWriterVTKBCMeshCases.cc \
# TestDataWriterPoints.cc \
# TestDataWriterVTKPoints.cc \
# TestDataWriterVTKPointsCases.cc \
# TestDataWriterFaultMesh.cc \
# TestDataWriterVTKFaultMesh.cc \
# TestDataWriterVTKFaultMeshCases.cc \
# TestOutputManager.cc \
# TestOutputSolnSubset.cc \
# TestOutputSolnPoints.cc \
# test_meshio.cc




noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh \
	TestOutputTriggerStep.hh \
	TestOutputTriggerTime.hh \
	FieldFactory.hh \
	TestOutputManager.hh \
	TestOutputSolnSubset.hh \
	TestOutputSolnPoints.hh \
	TestVertexFilterVecNorm.hh \
	TestDataWriter.hh \
	TestDataWriterMesh.hh \
	TestDataWriterVTK.hh \
	TestDataWriterVTKMesh.hh \
	TestDataWriterMaterial.hh \
	TestDataWriterVTKMaterial.hh \
	TestDataWriterSubmesh.hh \
	TestDataWriterVTKSubmesh.hh

# TestDataWriterSubmesh.hh \
# TestDataWriterBCMesh.hh \
# TestDataWriterFaultMesh.hh \
# TestDataWriterVTKMeshCases.hh \
# TestDataWriterVTKSubmesh.hh \
# TestDataWriterVTKSubmeshCases.hh \
# TestDataWriterVTKFaultMesh.hh \
# TestDataWriterVTKFaultMeshCases.hh \
# TestDataWriterVTKBCMesh.hh \
# TestDataWriterVTKBCMeshCases.hh \
# TestDataWriterPoints.hh \
# TestDataWriterVTKPoints.hh \
# TestDataWriterVTKPointsCases.hh



# Source files associated with testing data
# testmeshio_SOURCES += \
# 	data/DataWriterVTKDataSubmeshTri3.cc \
# 	data/DataWriterVTKDataSubmeshQuad4.cc \
# 	data/DataWriterVTKDataSubmeshTet4.cc \
# 	data/DataWriterVTKDataSubmeshHex8.cc \
# 	data/DataWriterVTKDataFaultMeshTri3.cc \
# 	data/DataWriterVTKDataFaultMeshQuad4.cc \
# 	data/DataWriterVTKDataFaultMeshTet4.cc \
# 	data/DataWriterVTKDataFaultMeshHex8.cc \
# 	data/DataWriterVTKDataBCMeshTri3.cc \
# 	data/DataWriterVTKDataBCMeshQuad4.cc \
# 	data/DataWriterVTKDataBCMeshTet4.cc \
# 	data/DataWriterVTKDataBCMeshHex8.cc \
# 	data/DataWriterHDF5DataSubmeshTri3.cc \
# 	data/DataWriterHDF5DataSubmeshQuad4.cc \
# 	data/DataWriterHDF5DataSubmeshTet4.cc \
# 	data/DataWriterHDF5DataSubmeshHex8.cc \
# 	data/DataWriterHDF5DataBCMeshTri3.cc \
# 	data/DataWriterHDF5DataBCMeshQuad4.cc \
# 	data/DataWriterHDF5DataBCMeshTet4.cc \
# 	data/DataWriterHDF5DataBCMeshHex8.cc \
# 	data/DataWriterHDF5DataFaultMeshTri3.cc \
# 	data/DataWriterHDF5DataFaultMeshQuad4.cc \
# 	data/DataWriterHDF5DataFaultMeshTet4.cc \
# 	data/DataWriterHDF5DataFaultMeshHex8.cc \
# 	data/OutputSolnPointsData.cc \
# 	data/OutputSolnPointsDataTri3.cc \
# 	data/OutputSolnPointsDataQuad4.cc \
# 	data/OutputSolnPointsDataTet4.cc \
# 	data/OutputSolnPointsDataHex8.cc \
# 	data/DataWriterDataPoints.cc \
# 	data/DataWriterVTKDataPointsTri3.cc \
# 	data/DataWriterVTKDataPointsQuad4.cc \
# 	data/DataWriterVTKDataPointsTet4.cc \
# 	data/DataWriterVTKDataPointsHex8.cc \
# 	data/DataWriterHDF5DataPointsTri3.cc \
# 	data/DataWriterHDF5DataPointsQuad4.cc \
# 	data/DataWriterHDF5DataPointsTet4.cc \
# 	data/DataWriterHDF5DataPointsHex8.cc

# noinst_HEADERS += \
# 	data/DataWriterVTKDataSubmeshTri3.hh \
# 	data/DataWriterVTKDataSubmeshQuad4.hh \
# 	data/DataWriterVTKDataSubmeshTet4.hh \
# 	data/DataWriterVTKDataSubmeshHex8.hh \
# 	data/DataWriterVTKDataFaultMeshTri3.hh \
# 	data/DataWriterVTKDataFaultMeshQuad4.hh \
# 	data/DataWriterVTKDataFaultMeshTet4.hh \
# 	data/DataWriterVTKDataFaultMeshHex8.hh \
# 	data/DataWriterVTKDataBCMeshTri3.hh \
# 	data/DataWriterVTKDataBCMeshQuad4.hh \
# 	data/DataWriterVTKDataBCMeshTet4.hh \
# 	data/DataWriterVTKDataBCMeshHex8.hh \
# 	data/DataWriterHDF5DataSubmeshTri3.hh \
# 	data/DataWriterHDF5DataSubmeshQuad4.hh \
# 	data/DataWriterHDF5DataSubmeshTet4.hh \
# 	data/DataWriterHDF5DataSubmeshHex8.hh \
# 	data/DataWriterHDF5DataBCMeshTri3.hh \
# 	data/DataWriterHDF5DataBCMeshQuad4.hh \
# 	data/DataWriterHDF5DataBCMeshTet4.hh \
# 	data/DataWriterHDF5DataBCMeshHex8.hh \
# 	data/DataWriterHDF5DataFaultMeshTri3.hh \
# 	data/DataWriterHDF5DataFaultMeshQuad4.hh \
# 	data/DataWriterHDF5DataFaultMeshTet4.hh \
# 	data/DataWriterHDF5DataFaultMeshHex8.hh \
# 	data/OutputSolnPointsData.hh \
# 	data/OutputSolnPointsDataTri3.hh \
# 	data/OutputSolnPointsDataQuad4.hh \
# 	data/OutputSolnPointsDataTet4.hh \
# 	data/OutputSolnPointsDataHex8.hh \
# 	data/DataWriterDataPoints.hh \
# 	data/DataWriterVTKDataPointsTri3.hh \
# 	data/DataWriterVTKDataPointsQuad4.hh \
# 	data/DataWriterVTKDataPointsTet4.hh \
# 	data/DataWriterVTKDataPointsHex8.hh \
# 	data/DataWriterHDF5DataPointsTri3.hh \
# 	data/DataWriterHDF5DataPointsQuad4.hh \
# 	data/DataWriterHDF5DataPointsTet4.hh \
# 	data/DataWriterHDF5DataPointsHex8.hh

AM_CPPFLAGS += \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testmeshio_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testmeshio_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testmeshio_SOURCES += \
	TestExodusII.cc \
	TestMeshIOCubit.cc \
	TestMeshIOCubit_Cases.cc
  noinst_HEADERS += \
	TestExodusII.hh \
	TestMeshIOCubit.hh
  testmeshio_LDADD += -lnetcdf
endif

if ENABLE_HDF5
  testmeshio_SOURCES += \
	TestHDF5.cc \
	TestDataWriterHDF5.cc \
	TestDataWriterHDF5Mesh.cc \
	TestDataWriterHDF5Mesh_Cases.cc \
	TestDataWriterHDF5ExtMesh.cc \
	TestDataWriterHDF5ExtMesh_Cases.cc 


#	TestDataWriterHDF5Material.cc \
#	TestDataWriterHDF5Material_Cases.cc \
#	TestDataWriterHDF5ExtMaterial.cc \
#	TestDataWriterHDF5ExtMaterial_Cases.cc


# TestDataWriterHDF5Submesh.cc \
# TestDataWriterHDF5SubmeshCases.cc \
# TestDataWriterHDF5Points.cc \
# TestDataWriterHDF5PointsCases.cc \
# TestDataWriterHDF5BCMesh.cc \
# TestDataWriterHDF5BCMeshCases.cc \
# TestDataWriterHDF5FaultMesh.cc \
# TestDataWriterHDF5FaultMeshCases.cc \
# TestDataWriterHDF5ExtSubmesh.cc \
# TestDataWriterHDF5ExtSubmeshCases.cc \
# TestDataWriterHDF5ExtPoints.cc \
# TestDataWriterHDF5ExtPointsCases.cc \
# TestDataWriterHDF5ExtBCMesh.cc \
# TestDataWriterHDF5ExtBCMeshCases.cc \
# TestDataWriterHDF5ExtFaultMesh.cc \
# TestDataWriterHDF5ExtFaultMeshCases.cc


  noinst_HEADERS += \
	TestHDF5.hh \
	TestDataWriterHDF5.hh \
	TestDataWriterHDF5Mesh.hh \
	TestDataWriterHDF5Material.hh \
	TestDataWriterHDF5ExtMesh.hh \
	TestDataWriterHDF5ExtMaterial.hh

# TestDataWriterHDF5Submesh.hh \
# TestDataWriterHDF5SubmeshCases.hh \
# TestDataWriterHDF5Points.hh \
# TestDataWriterHDF5PointsCases.hh \
# TestDataWriterHDF5BCMesh.hh \
# TestDataWriterHDF5BCMeshCases.hh \
# TestDataWriterHDF5FaultMesh.hh \
# TestDataWriterHDF5FaultMeshCases.hh \
# TestDataWriterHDF5ExtMesh.hh \
# TestDataWriterHDF5ExtSubmesh.hh \
# TestDataWriterHDF5ExtSubmeshCases.hh \
# TestDataWriterHDF5ExtPoints.hh \
# TestDataWriterHDF5ExtPointsCases.hh \
# TestDataWriterHDF5ExtBCMesh.hh \
# TestDataWriterHDF5ExtBCMeshCases.hh \
# TestDataWriterHDF5ExtFaultMesh.hh \
# TestDataWriterHDF5ExtFaultMeshCases.hh


  testmeshio_LDADD += -lhdf5
endif


noinst_tmp = \
	mesh1D.txt \
	mesh1Din2D.txt \
	mesh1Din3D.txt \
	mesh2D.txt \
	mesh2Din3D.txt \
	mesh3D.txt


CLEANFILES = $(noinst_tmp)

clean-local: clean-local-tmp
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.vtk *.dat *.dat.info *.h5 *.xmf

# End of file
