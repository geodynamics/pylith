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
	TestMeshIO.cc \
	TestMeshIOAscii.cc \
	TestMeshIOLagrit.cc \
	test_meshio.cc

noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh

# Source files associated with testing data
testmeshio_SOURCES += \
	data/MeshData.cc \
	data/MeshData1D.cc \
	data/MeshData1Din3D.cc \
	data/MeshData1Din2D.cc \
	data/MeshData2D.cc \
	data/MeshData2Din3D.cc \
	data/MeshData3D.cc \
	data/MeshDataCubitTri.cc \
	data/MeshDataCubitQuad.cc \
	data/MeshDataCubitTet.cc \
	data/MeshDataCubitHex.cc \
	data/MeshDataLagritTet.cc

noinst_HEADERS += \
	data/MeshData.hh \
	data/MeshData1D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData2D.hh \
	data/MeshData2Din3D.hh \
	data/MeshData3D.hh \
	data/MeshDataCubitTri.hh \
	data/MeshDataCubitQuad.hh \
	data/MeshDataCubitTet.hh \
	data/MeshDataCubitHex.hh \
	data/MeshDataLagritTet.hh

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

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
	mesh3D.txt

CLEANFILES = $(noinst_tmp)

# End of file 
