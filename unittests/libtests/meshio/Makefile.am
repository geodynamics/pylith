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
# Copyright (c) 2010-2014 University of California, Davis
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
	TestMeshIOLagrit.cc \
	TestCellFilterAvg.cc \
	TestVertexFilterVecNorm.cc \
	TestDataWriterMesh.cc \
	TestDataWriterSubMesh.cc \
	TestDataWriterBCMesh.cc \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMeshCases.cc \
	TestDataWriterVTKMatMeshCases.cc \
	TestDataWriterVTKSubMesh.cc \
	TestDataWriterVTKSubMeshCases.cc \
	TestDataWriterVTKBCMesh.cc \
	TestDataWriterVTKBCMeshCases.cc \
	TestDataWriterPoints.cc \
	TestDataWriterVTKPoints.cc \
	TestDataWriterVTKPointsCases.cc \
	TestDataWriterFaultMesh.cc \
	TestDataWriterVTKFaultMesh.cc \
	TestDataWriterVTKFaultMeshCases.cc \
	TestOutputManager.cc \
	TestOutputSolnSubset.cc \
	TestOutputSolnPoints.cc \
	test_meshio.cc



noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh \
	TestOutputManager.hh \
	TestOutputSolnSubset.hh \
	TestOutputSolnPoints.hh \
	TestVertexFilterVecNorm.hh \
	TestCellFilterAvg.hh \
	TestDataWriterMesh.hh \
	TestDataWriterVTK.hh \
	TestDataWriterSubMesh.hh \
	TestDataWriterBCMesh.hh \
	TestDataWriterFaultMesh.hh \
	TestDataWriterVTKMesh.hh \
	TestDataWriterVTKMeshCases.hh \
	TestDataWriterVTKMatMeshCases.hh \
	TestDataWriterVTKSubMesh.hh \
	TestDataWriterVTKSubMeshCases.hh \
	TestDataWriterVTKFaultMesh.hh \
	TestDataWriterVTKFaultMeshCases.hh \
	TestDataWriterVTKBCMesh.hh \
	TestDataWriterVTKBCMeshCases.hh \
	TestDataWriterPoints.hh \
	TestDataWriterVTKPoints.hh \
	TestDataWriterVTKPointsCases.hh


# Source files associated with testing data
testmeshio_SOURCES += \
	data/DataWriterData.cc \
	data/DataWriterVTKDataMeshTri3.cc \
	data/DataWriterVTKDataMeshQuad4.cc \
	data/DataWriterVTKDataMeshTet4.cc \
	data/DataWriterVTKDataMeshHex8.cc \
	data/DataWriterVTKDataMatMeshTri3.cc \
	data/DataWriterVTKDataMatMeshQuad4.cc \
	data/DataWriterVTKDataMatMeshTet4.cc \
	data/DataWriterVTKDataMatMeshHex8.cc \
	data/DataWriterVTKDataSubMeshTri3.cc \
	data/DataWriterVTKDataSubMeshQuad4.cc \
	data/DataWriterVTKDataSubMeshTet4.cc \
	data/DataWriterVTKDataSubMeshHex8.cc \
	data/DataWriterVTKDataFaultMeshTri3.cc \
	data/DataWriterVTKDataFaultMeshQuad4.cc \
	data/DataWriterVTKDataFaultMeshTet4.cc \
	data/DataWriterVTKDataFaultMeshHex8.cc \
	data/DataWriterVTKDataBCMeshTri3.cc \
	data/DataWriterVTKDataBCMeshQuad4.cc \
	data/DataWriterVTKDataBCMeshTet4.cc \
	data/DataWriterVTKDataBCMeshHex8.cc \
	data/DataWriterHDF5DataMeshTri3.cc \
	data/DataWriterHDF5DataMeshQuad4.cc \
	data/DataWriterHDF5DataMeshTet4.cc \
	data/DataWriterHDF5DataMeshHex8.cc \
	data/DataWriterHDF5DataMatMeshTri3.cc \
	data/DataWriterHDF5DataMatMeshQuad4.cc \
	data/DataWriterHDF5DataMatMeshTet4.cc \
	data/DataWriterHDF5DataMatMeshHex8.cc \
	data/DataWriterHDF5DataSubMeshTri3.cc \
	data/DataWriterHDF5DataSubMeshQuad4.cc \
	data/DataWriterHDF5DataSubMeshTet4.cc \
	data/DataWriterHDF5DataSubMeshHex8.cc \
	data/DataWriterHDF5DataBCMeshTri3.cc \
	data/DataWriterHDF5DataBCMeshQuad4.cc \
	data/DataWriterHDF5DataBCMeshTet4.cc \
	data/DataWriterHDF5DataBCMeshHex8.cc \
	data/DataWriterHDF5DataFaultMeshTri3.cc \
	data/DataWriterHDF5DataFaultMeshQuad4.cc \
	data/DataWriterHDF5DataFaultMeshTet4.cc \
	data/DataWriterHDF5DataFaultMeshHex8.cc \
	data/OutputSolnPointsData.cc \
	data/OutputSolnPointsDataTri3.cc \
	data/OutputSolnPointsDataQuad4.cc \
	data/OutputSolnPointsDataTet4.cc \
	data/OutputSolnPointsDataHex8.cc \
	data/DataWriterDataPoints.cc \
	data/DataWriterVTKDataPointsTri3.cc \
	data/DataWriterVTKDataPointsQuad4.cc \
	data/DataWriterVTKDataPointsTet4.cc \
	data/DataWriterVTKDataPointsHex8.cc \
	data/DataWriterHDF5DataPointsTri3.cc \
	data/DataWriterHDF5DataPointsQuad4.cc \
	data/DataWriterHDF5DataPointsTet4.cc \
	data/DataWriterHDF5DataPointsHex8.cc \
	data/MeshData.cc \
	data/MeshData1D.cc \
	data/MeshData1Din3D.cc \
	data/MeshData1Din2D.cc \
	data/MeshData2D.cc \
	data/MeshData2Din3D.cc \
	data/MeshData3D.cc \
	data/MeshData3DIndexOne.cc \
	data/MeshDataCubitTri.cc \
	data/MeshDataCubitQuad.cc \
	data/MeshDataCubitTet.cc \
	data/MeshDataCubitHex.cc \
	data/MeshDataLagritTet.cc


noinst_HEADERS += \
	data/DataWriterData.hh \
	data/DataWriterVTKDataMeshTri3.hh \
	data/DataWriterVTKDataMeshQuad4.hh \
	data/DataWriterVTKDataMeshTet4.hh \
	data/DataWriterVTKDataMeshHex8.hh \
	data/DataWriterVTKDataMatMeshTri3.hh \
	data/DataWriterVTKDataMatMeshQuad4.hh \
	data/DataWriterVTKDataMatMeshTet4.hh \
	data/DataWriterVTKDataMatMeshHex8.hh \
	data/DataWriterVTKDataSubMeshTri3.hh \
	data/DataWriterVTKDataSubMeshQuad4.hh \
	data/DataWriterVTKDataSubMeshTet4.hh \
	data/DataWriterVTKDataSubMeshHex8.hh \
	data/DataWriterVTKDataFaultMeshTri3.hh \
	data/DataWriterVTKDataFaultMeshQuad4.hh \
	data/DataWriterVTKDataFaultMeshTet4.hh \
	data/DataWriterVTKDataFaultMeshHex8.hh \
	data/DataWriterVTKDataBCMeshTri3.hh \
	data/DataWriterVTKDataBCMeshQuad4.hh \
	data/DataWriterVTKDataBCMeshTet4.hh \
	data/DataWriterVTKDataBCMeshHex8.hh \
	data/DataWriterHDF5DataMeshTri3.hh \
	data/DataWriterHDF5DataMeshQuad4.hh \
	data/DataWriterHDF5DataMeshTet4.hh \
	data/DataWriterHDF5DataMeshHex8.hh \
	data/DataWriterHDF5DataMatMeshTri3.hh \
	data/DataWriterHDF5DataMatMeshQuad4.hh \
	data/DataWriterHDF5DataMatMeshTet4.hh \
	data/DataWriterHDF5DataMatMeshHex8.hh \
	data/DataWriterHDF5DataSubMeshTri3.hh \
	data/DataWriterHDF5DataSubMeshQuad4.hh \
	data/DataWriterHDF5DataSubMeshTet4.hh \
	data/DataWriterHDF5DataSubMeshHex8.hh \
	data/DataWriterHDF5DataBCMeshTri3.hh \
	data/DataWriterHDF5DataBCMeshQuad4.hh \
	data/DataWriterHDF5DataBCMeshTet4.hh \
	data/DataWriterHDF5DataBCMeshHex8.hh \
	data/DataWriterHDF5DataFaultMeshTri3.hh \
	data/DataWriterHDF5DataFaultMeshQuad4.hh \
	data/DataWriterHDF5DataFaultMeshTet4.hh \
	data/DataWriterHDF5DataFaultMeshHex8.hh \
	data/OutputSolnPointsData.hh \
	data/OutputSolnPointsDataTri3.hh \
	data/OutputSolnPointsDataQuad4.hh \
	data/OutputSolnPointsDataTet4.hh \
	data/OutputSolnPointsDataHex8.hh \
	data/DataWriterDataPoints.hh \
	data/DataWriterVTKDataPointsTri3.hh \
	data/DataWriterVTKDataPointsQuad4.hh \
	data/DataWriterVTKDataPointsTet4.hh \
	data/DataWriterVTKDataPointsHex8.hh \
	data/DataWriterHDF5DataPointsTri3.hh \
	data/DataWriterHDF5DataPointsQuad4.hh \
	data/DataWriterHDF5DataPointsTet4.hh \
	data/DataWriterHDF5DataPointsHex8.hh \
	data/MeshData.hh \
	data/MeshData1D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData1Din3D.hh \
	data/MeshData2D.hh \
	data/MeshData2Din3D.hh \
	data/MeshData3D.hh \
	data/MeshData3DIndexOne.hh \
	data/MeshDataCubitTri.hh \
	data/MeshDataCubitQuad.hh \
	data/MeshDataCubitTet.hh \
	data/MeshDataCubitHex.hh \
	data/MeshDataLagritTet.hh


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
	TestMeshIOCubit.cc
  noinst_HEADERS += \
	TestExodusII.hh \
	TestMeshIOCubit.hh
  testmeshio_LDADD += -lnetcdf_c++ -lnetcdf
endif

if ENABLE_HDF5
  testmeshio_SOURCES += \
	TestHDF5.cc \
	TestDataWriterHDF5.cc \
	TestDataWriterHDF5Mesh.cc \
	TestDataWriterHDF5MeshCases.cc \
	TestDataWriterHDF5MatMeshCases.cc \
	TestDataWriterHDF5SubMesh.cc \
	TestDataWriterHDF5SubMeshCases.cc \
	TestDataWriterHDF5Points.cc \
	TestDataWriterHDF5PointsCases.cc \
	TestDataWriterHDF5BCMesh.cc \
	TestDataWriterHDF5BCMeshCases.cc \
	TestDataWriterHDF5FaultMesh.cc \
	TestDataWriterHDF5FaultMeshCases.cc \
	TestDataWriterHDF5ExtMesh.cc \
	TestDataWriterHDF5ExtMeshCases.cc \
	TestDataWriterHDF5ExtMatMeshCases.cc \
	TestDataWriterHDF5ExtSubMesh.cc \
	TestDataWriterHDF5ExtSubMeshCases.cc \
	TestDataWriterHDF5ExtPoints.cc \
	TestDataWriterHDF5ExtPointsCases.cc \
	TestDataWriterHDF5ExtBCMesh.cc \
	TestDataWriterHDF5ExtBCMeshCases.cc \
	TestDataWriterHDF5ExtFaultMesh.cc \
	TestDataWriterHDF5ExtFaultMeshCases.cc \
	TestXdmf.cc

  noinst_HEADERS += \
	TestHDF5.hh \
	TestDataWriterHDF5.hh \
	TestDataWriterHDF5Mesh.hh \
	TestDataWriterHDF5MeshCases.hh \
	TestDataWriterHDF5MatMeshCases.hh \
	TestDataWriterHDF5SubMesh.hh \
	TestDataWriterHDF5SubMeshCases.hh \
	TestDataWriterHDF5Points.hh \
	TestDataWriterHDF5PointsCases.hh \
	TestDataWriterHDF5BCMesh.hh \
	TestDataWriterHDF5BCMeshCases.hh \
	TestDataWriterHDF5FaultMesh.hh \
	TestDataWriterHDF5FaultMeshCases.hh \
	TestDataWriterHDF5ExtMesh.hh \
	TestDataWriterHDF5ExtMeshCases.hh \
	TestDataWriterHDF5ExtMatMeshCases.hh \
	TestDataWriterHDF5ExtSubMesh.hh \
	TestDataWriterHDF5ExtSubMeshCases.hh \
	TestDataWriterHDF5ExtPoints.hh \
	TestDataWriterHDF5ExtPointsCases.hh \
	TestDataWriterHDF5ExtBCMesh.hh \
	TestDataWriterHDF5ExtBCMeshCases.hh \
	TestDataWriterHDF5ExtFaultMesh.hh \
	TestDataWriterHDF5ExtFaultMeshCases.hh \
	TestXdmf.hh

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


leakcheck: testmeshio
	valgrind --log-file=valgrind_meshio.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testmeshio


# End of file 
