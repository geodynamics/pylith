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

SUBDIRS = data

TESTS = testmeshio

check_PROGRAMS = testmeshio

# Primary source files
testmeshio_SOURCES = \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMeshLine2.cc \
	TestDataWriterVTKMeshTri3.cc \
	TestDataWriterVTKMeshQuad4.cc \
	TestDataWriterVTKMeshTet4.cc \
	TestDataWriterVTKMeshHex8.cc \
	TestDataWriterVTKSubMeshLine2.cc \
	TestDataWriterVTKSubMeshTri3.cc \
	TestDataWriterVTKSubMeshQuad4.cc \
	TestDataWriterVTKSubMeshTet4.cc \
	TestDataWriterVTKSubMeshHex8.cc \
	TestMeshIO.cc \
	TestMeshIOAscii.cc \
	TestMeshIOLagrit.cc \
	test_meshio.cc

noinst_HEADERS = \
	TestDataWriterVTK.hh \
	TestDataWriterVTKMesh.hh \
	TestDataWriterVTKMeshLine2.hh \
	TestDataWriterVTKMeshTri3.hh \
	TestDataWriterVTKMeshQuad4.hh \
	TestDataWriterVTKMeshTet4.hh \
	TestDataWriterVTKMeshHex8.hh \
	TestDataWriterVTKSubMeshLine2.hh \
	TestDataWriterVTKSubMeshTri3.hh \
	TestDataWriterVTKSubMeshQuad4.hh \
	TestDataWriterVTKSubMeshTet4.hh \
	TestDataWriterVTKSubMeshHex8.hh \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh

# Source files associated with testing data
testmeshio_SOURCES += \
	data/DataWriterVTKData.cc \
	data/DataWriterVTKDataMesh.cc \
	data/DataWriterVTKDataMeshLine2.cc \
	data/DataWriterVTKDataMeshTri3.cc \
	data/DataWriterVTKDataMeshQuad4.cc \
	data/DataWriterVTKDataMeshTet4.cc \
	data/DataWriterVTKDataMeshHex8.cc \
	data/DataWriterVTKDataSubMeshLine2.cc \
	data/DataWriterVTKDataSubMeshTri3.cc \
	data/DataWriterVTKDataSubMeshQuad4.cc \
	data/DataWriterVTKDataSubMeshTet4.cc \
	data/DataWriterVTKDataSubMeshHex8.cc \
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
	data/DataWriterVTKData.hh \
	data/DataWriterVTKDataMesh.hh \
	data/DataWriterVTKDataMeshLine2.hh \
	data/DataWriterVTKDataMeshTri3.hh \
	data/DataWriterVTKDataMeshQuad4.hh \
	data/DataWriterVTKDataMeshTet4.hh \
	data/DataWriterVTKDataMeshHex8.hh \
	data/DataWriterVTKDataSubMeshLine2.hh \
	data/DataWriterVTKDataSubMeshTri3.hh \
	data/DataWriterVTKDataSubMeshQuad4.hh \
	data/DataWriterVTKDataSubMeshTet4.hh \
	data/DataWriterVTKDataSubMeshHex8.hh \
	data/MeshData.hh \
	data/MeshData1D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData1Din2D.hh \
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
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
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
	line2_sub_t10.vtk \
	line2_sub_vertex_t10.vtk \
	line2_sub_cell_t10.vtk \
	tri3_sub_t10.vtk \
	tri3_sub_vertex_t10.vtk \
	tri3_sub_cell_t10.vtk \
	quad4_sub_t10.vtk \
	quad4_sub_vertex_t10.vtk \
	quad4_sub_cell_t10.vtk \
	tet4_sub_t10.vtk \
	tet4_sub_vertex_t10.vtk \
	tet4_sub_cell_t10.vtk \
	hex8_sub_t10.vtk \
	hex8_sub_vertex_t10.vtk \
	hex8_sub_cell_t10.vtk

CLEANFILES = $(noinst_tmp)

# End of file 
