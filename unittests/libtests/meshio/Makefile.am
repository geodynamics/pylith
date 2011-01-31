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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = meshio
include $(top_srcdir)/subpackage.am

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
	TestDataWriterFaultMesh.cc \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMeshCases.cc \
	TestDataWriterVTKMatMeshCases.cc \
	TestDataWriterVTKBCMesh.cc \
	TestDataWriterVTKBCMeshCases.cc \
	TestDataWriterVTKSubMesh.cc \
	TestDataWriterVTKSubMeshCases.cc \
	TestOutputManager.cc \
	TestOutputSolnSubset.cc \
	TestDataWriterVTKFaultMesh.cc \
	TestDataWriterVTKFaultMeshCases.cc \
	test_meshio.cc


noinst_HEADERS = \
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
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh \
	TestOutputManager.hh \
	TestOutputSolnSubset.hh \
	TestVertexFilterVecNorm.hh


# Source files associated with testing data
testmeshio_SOURCES += \
	data/DataWriterData.cc \
	data/DataWriterVTKDataMeshLine2.cc \
	data/DataWriterVTKDataMeshTri3.cc \
	data/DataWriterVTKDataMeshQuad4.cc \
	data/DataWriterVTKDataMeshTet4.cc \
	data/DataWriterVTKDataMeshHex8.cc \
	data/DataWriterVTKDataMatMeshLine2.cc \
	data/DataWriterVTKDataMatMeshTri3.cc \
	data/DataWriterVTKDataMatMeshQuad4.cc \
	data/DataWriterVTKDataMatMeshTet4.cc \
	data/DataWriterVTKDataMatMeshHex8.cc \
	data/DataWriterVTKDataSubMeshLine2.cc \
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
	data/DataWriterHDF5DataMeshLine2.cc \
	data/DataWriterHDF5DataMeshTri3.cc \
	data/DataWriterHDF5DataMeshQuad4.cc \
	data/DataWriterHDF5DataMeshTet4.cc \
	data/DataWriterHDF5DataMeshHex8.cc \
	data/DataWriterHDF5DataMatMeshLine2.cc \
	data/DataWriterHDF5DataMatMeshTri3.cc \
	data/DataWriterHDF5DataMatMeshQuad4.cc \
	data/DataWriterHDF5DataMatMeshTet4.cc \
	data/DataWriterHDF5DataMatMeshHex8.cc \
	data/DataWriterHDF5DataSubMeshLine2.cc \
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
	data/DataWriterVTKDataMeshLine2.hh \
	data/DataWriterVTKDataMeshTri3.hh \
	data/DataWriterVTKDataMeshQuad4.hh \
	data/DataWriterVTKDataMeshTet4.hh \
	data/DataWriterVTKDataMeshHex8.hh \
	data/DataWriterVTKDataMatMeshLine2.hh \
	data/DataWriterVTKDataMatMeshTri3.hh \
	data/DataWriterVTKDataMatMeshQuad4.hh \
	data/DataWriterVTKDataMatMeshTet4.hh \
	data/DataWriterVTKDataMatMeshHex8.hh \
	data/DataWriterVTKDataSubMeshLine2.hh \
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
	data/DataWriterHDF5DataMeshLine2.hh \
	data/DataWriterHDF5DataMeshTri3.hh \
	data/DataWriterHDF5DataMeshQuad4.hh \
	data/DataWriterHDF5DataMeshTet4.hh \
	data/DataWriterHDF5DataMeshHex8.hh \
	data/DataWriterHDF5DataMatMeshLine2.hh \
	data/DataWriterHDF5DataMatMeshTri3.hh \
	data/DataWriterHDF5DataMatMeshQuad4.hh \
	data/DataWriterHDF5DataMatMeshTet4.hh \
	data/DataWriterHDF5DataMatMeshHex8.hh \
	data/DataWriterHDF5DataSubMeshLine2.hh \
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


AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testmeshio_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testmeshio_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testmeshio_SOURCES += TestMeshIOCubit.cc
  noinst_HEADERS += TestMeshIOCubit.hh
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
	TestDataWriterHDF5BCMesh.cc \
	TestDataWriterHDF5BCMeshCases.cc \
	TestDataWriterHDF5FaultMesh.cc \
	TestDataWriterHDF5FaultMeshCases.cc \
	TestDataWriterHDF5ExtMesh.cc \
	TestDataWriterHDF5ExtMeshCases.cc \
	TestDataWriterHDF5ExtMatMeshCases.cc

  noinst_HEADERS += \
	TestDataWriterHDF5.hh \
	TestDataWriterHDF5Mesh.hh \
	TestDataWriterHDF5MeshCases.hh \
	TestDataWriterHDF5MatMeshCases.hh \
	TestDataWriterHDF5SubMesh.hh \
	TestDataWriterHDF5SubMeshCases.hh \
	TestDataWriterHDF5BCMesh.hh \
	TestDataWriterHDF5BCMeshCases.hh \
	TestDataWriterHDF5FaultMesh.hh \
	TestDataWriterHDF5FaultMeshCases.hh \
	TestDataWriterHDF5ExtMesh.hh \
	TestDataWriterHDF5ExtMeshCases.hh \
	TestDataWriterHDF5ExtMatMeshCases.hh

  testmeshio_LDADD += -lhdf5
endif


noinst_tmp = \
	mesh1D.txt \
	mesh1Din2D.txt \
	mesh1Din3D.txt \
	mesh2D.txt \
	mesh2Din3D.txt \
	mesh3D.txt \
	line2_t10.vtk \
	line2_vertex_t10.vtk \
	line2_cell_t10.vtk \
	tri3_t10.vtk \
	tri3_vertex_t10.vtk \
	tri3_cell_t10.vtk \
	quad4_t10.vtk \
	quad4_vertex_t10.vtk \
	quad4_cell_t10.vtk \
	tet4_t10.vtk \
	tet4_vertex_t10.vtk \
	tet4_cell_t10.vtk \
	hex8_t10.vtk \
	hex8_vertex_t10.vtk \
	hex8_cell_t10.vtk \
	line2_mat_t10.vtk \
	line2_mat_vertex_t10.vtk \
	line2_mat_cell_t10.vtk \
	tri3_mat_t10.vtk \
	tri3_mat_vertex_t10.vtk \
	tri3_mat_cell_t10.vtk \
	quad4_mat_t10.vtk \
	quad4_mat_vertex_t10.vtk \
	quad4_mat_cell_t10.vtk \
	tet4_mat_t10.vtk \
	tet4_mat_vertex_t10.vtk \
	tet4_mat_cell_t10.vtk \
	hex8_mat_t10.vtk \
	hex8_mat_vertex_t10.vtk \
	hex8_mat_cell_t10.vtk \
	line2_surf_t10.vtk \
	line2_surf_vertex_t10.vtk \
	line2_surf_cell_t10.vtk \
	tri3_surf_t10.vtk \
	tri3_surf_vertex_t10.vtk \
	tri3_surf_cell_t10.vtk \
	quad4_surf_t10.vtk \
	quad4_surf_vertex_t10.vtk \
	quad4_surf_cell_t10.vtk \
	tet4_surf_t10.vtk \
	tet4_surf_vertex_t10.vtk \
	tet4_surf_cell_t10.vtk \
	hex8_surf_t10.vtk \
	hex8_surf_vertex_t10.vtk \
	hex8_surf_cell_t10.vtk \
	tri3_fault_t10.vtk \
	tri3_fault_vertex_t10.vtk \
	tri3_fault_cell_t10.vtk \
	quad4_fault_t10.vtk \
	quad4_fault_vertex_t10.vtk \
	quad4_fault_cell_t10.vtk \
	tet4_fault_t10.vtk \
	tet4_fault_vertex_t10.vtk \
	tet4_fault_cell_t10.vtk \
	hex8_fault_t10.vtk \
	hex8_fault_vertex_t10.vtk \
	hex8_fault_cell_t10.vtk \
	tri3_bc_t10.vtk \
	tri3_bc_vertex_t10.vtk \
	tri3_bc_cell_t10.vtk \
	quad4_bc_t10.vtk \
	quad4_bc_vertex_t10.vtk \
	quad4_bc_cell_t10.vtk \
	tet4_bc_t10.vtk \
	tet4_bc_vertex_t10.vtk \
	tet4_bc_cell_t10.vtk \
	hex8_bc_t10.vtk \
	hex8_bc_vertex_t10.vtk \
	hex8_bc_cell_t10.vtk \
	output_t12.vtk \
	output_vertex_t12.vtk \
	output_vertex_filter_t12.vtk \
	output_cell_t12.vtk \
	output_cell_filter_t12.vtk \
	test.h5 \
	test.dat \
	hex8.h5 \
	hex8_vertex.h5 \
	hex8_cell.h5

CLEANFILES = $(noinst_tmp)

# End of file 
