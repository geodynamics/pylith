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
	TestMeshIOCubit.cc \
	test_meshio.cc

noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOCubit.hh

# Source files associated with testing data
testmeshio_SOURCES += \
	data/MeshData.cc \
	data/MeshData1D.cc \
	data/MeshData1Din3D.cc \
	data/MeshData1Din2D.cc \
	data/MeshData2D.cc \
	data/MeshData2Din3D.cc \
	data/MeshData3D.cc \
	data/MeshDataCubitHex.cc \
	data/MeshDataCubitTet.cc

noinst_HEADERS += \
	data/MeshData.hh \
	data/MeshData1D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData1Din2D.hh \
	data/MeshData2D.hh \
	data/MeshData2Din3D.hh \
	data/MeshData3D.hh \
	data/MeshDataCubitHex.hh \
	data/MeshDataCubitTet.hh

testmeshio_LDFLAGS = $(PETSC_LIB)
#$(PYTHON_BLDLIBRARY)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testmeshio_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la

if ENABLE_CUBIT
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
