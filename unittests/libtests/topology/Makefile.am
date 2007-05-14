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

subpackage = topology
include $(top_srcdir)/subpackage.am

TESTS = testtopology

check_PROGRAMS = testtopology

# Primary source files
testtopology_SOURCES = \
	TestCellGeometry.cc \
	TestGeometryPoint1D.cc \
	TestGeometryPoint2D.cc \
	TestGeometryPoint3D.cc \
	TestGeometryLine1D.cc \
	TestGeometryLine2D.cc \
	TestGeometryLine3D.cc \
	TestGeometryTri2D.cc \
	TestGeometryTri3D.cc \
	TestGeometryQuad2D.cc \
	TestGeometryTet3D.cc \
	test_topology.cc

noinst_HEADERS =
	TestCellGeometry.hh \
	TestGeometryPoint1D.hh \
	TestGeometryPoint2D.hh \
	TestGeometryPoint3D.hh \
	TestGeometryLine1D.hh \
	TestGeometryLine2D.hh \
	TestGeometryLine3D.hh \
	TestGeometryTri2D.hh \
	TestGeometryTri3D.hh \
	TestGeometryQuad2D.hh \
	TestGeometryTet3D.hh

# Source files associated with testing data
testtopology_SOURCES += \
	data/CellGeomData.cc \
	data/GeomDataPoint1D.cc \
	data/GeomDataPoint2D.cc \
	data/GeomDataPoint3D.cc \
	data/GeomDataLine1D.cc \
	data/GeomDataLine2D.cc \
	data/GeomDataLine3D.cc \
	data/GeomDataTri2D.cc \
	data/GeomDataTri3D.cc \
	data/GeomDataQuad2D.cc \
	data/GeomDataTet3D.cc

noinst_HEADERS += 
	data/CellGeomData.hh \
	data/GeometryDataPoint1D.hh \
	data/GeometryDataPoint2D.hh \
	data/GeometryDataPoint3D.hh \
	data/GeometryDataLine1D.hh \
	data/GeometryDataLine2D.hh \
	data/GeometryDataLine3D.hh \
	data/GeometryDataTri2D.hh \
	data/GeometryDataTri3D.hh \
	data/GeometryDataQuad2D.hh \
	data/GeometryDataTet3D.hh

testtopology_LDFLAGS = $(PETSC_LIB) $(PYTHON_BLDLIBRARY)

AM_CPPFLAGS = \
	$(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testtopology_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la

if ENABLE_CUBIT
  testtopology_LDADD += -lnetcdf_c++ -lnetcdf
endif

noinst_tmp =

CLEANFILES = $(noinst_tmp)

# End of file 
