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

subpackage = feassemble
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testfeassemble

check_PROGRAMS = testfeassemble

# Primary source files
testfeassemble_SOURCES = \
	TestCellGeometry.cc \
	TestGeometryLine2D.cc \
	TestGeometryLine3D.cc \
	TestGeometryTri2D.cc \
	TestGeometryTri3D.cc \
	TestGeometryQuad2D.cc \
	TestGeometryQuad3D.cc \
	TestGeometryTet3D.cc \
	TestGeometryHex3D.cc \
	TestQuadratureRefCell.cc \
	TestQuadratureEngine.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature1Din3D.cc \
	TestQuadrature2D.cc \
	TestQuadrature2Din3D.cc \
	TestQuadrature3D.cc \
	TestQuadrature.cc \
	TestIntegrator.cc \
	TestIntegratorElasticity.cc \
	TestElasticityExplicit.cc \
	TestElasticityExplicitCases.cc \
	TestElasticityExplicitTri3.cc \
	TestElasticityExplicitTet4.cc \
	TestElasticityImplicit.cc \
	TestElasticityImplicitCases.cc \
	TestIntegratorElasticityLgDeform.cc \
	TestElasticityExplicitLgDeform.cc \
	TestElasticityExplicitLgDeformCases.cc \
	TestElasticityImplicitLgDeform.cc \
	TestElasticityImplicitLgDeformCases.cc \
	test_feassemble.cc



noinst_HEADERS = \
	TestCellGeometry.hh \
	TestElasticityExplicit.hh \
	TestElasticityExplicitCases.hh \
	TestElasticityExplicitTri3.hh \
	TestElasticityExplicitTet4.hh \
	TestElasticityImplicit.hh \
	TestElasticityImplicitCases.hh \
	TestIntegratorElasticityLgDeform.hh \
	TestElasticityExplicitLgDeform.hh \
	TestElasticityExplicitLgDeformCases.hh \
	TestElasticityImplicitLgDeform.hh \
	TestElasticityImplicitLgDeformCases.hh \
	TestGeometryLine2D.hh \
	TestGeometryLine3D.hh \
	TestGeometryTri2D.hh \
	TestGeometryTri3D.hh \
	TestGeometryQuad2D.hh \
	TestGeometryQuad3D.hh \
	TestGeometryTet3D.hh \
	TestGeometryHex3D.hh \
	TestQuadratureRefCell.hh \
	TestQuadratureEngine.hh \
	TestIntegrator.hh \
	TestIntegratorElasticity.hh \
	TestQuadrature.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature1Din3D.hh \
	TestQuadrature2D.hh \
	TestQuadrature2Din3D.hh \
	TestQuadrature3D.hh

# Source files associated with testing data
testfeassemble_SOURCES += \
	data/CellGeomData.cc \
	data/GeomDataLine2D.cc \
	data/GeomDataLine3D.cc \
	data/GeomDataTri2D.cc \
	data/GeomDataTri3D.cc \
	data/GeomDataQuad2D.cc \
	data/GeomDataQuad3D.cc \
	data/GeomDataTet3D.cc \
	data/GeomDataHex3D.cc \
	data/IntegratorData.cc \
	data/ElasticityExplicitData.cc \
	data/ElasticityExplicitData2DLinear.cc \
	data/ElasticityExplicitData2DQuadratic.cc \
	data/ElasticityExplicitData3DLinear.cc \
	data/ElasticityExplicitData3DQuadratic.cc \
	data/ElasticityExplicitGravData2DLinear.cc \
	data/ElasticityExplicitGravData2DQuadratic.cc \
	data/ElasticityExplicitGravData3DLinear.cc \
	data/ElasticityExplicitGravData3DQuadratic.cc \
	data/ElasticityImplicitData2DLinear.cc \
	data/ElasticityImplicitData2DQuadratic.cc \
	data/ElasticityImplicitData3DLinear.cc \
	data/ElasticityImplicitData3DQuadratic.cc \
	data/ElasticityImplicitGravData2DLinear.cc \
	data/ElasticityImplicitGravData2DQuadratic.cc \
	data/ElasticityImplicitGravData3DLinear.cc \
	data/ElasticityImplicitGravData3DQuadratic.cc \
	data/ElasticityExplicitLgDeformData2DLinear.cc \
	data/ElasticityExplicitLgDeformData2DQuadratic.cc \
	data/ElasticityExplicitLgDeformData3DLinear.cc \
	data/ElasticityExplicitLgDeformData3DQuadratic.cc \
	data/ElasticityExplicitLgDeformGravData2DLinear.cc \
	data/ElasticityExplicitLgDeformGravData2DQuadratic.cc \
	data/ElasticityExplicitLgDeformGravData3DLinear.cc \
	data/ElasticityExplicitLgDeformGravData3DQuadratic.cc \
	data/ElasticityImplicitLgDeformData2DLinear.cc \
	data/ElasticityImplicitLgDeformData2DQuadratic.cc \
	data/ElasticityImplicitLgDeformData3DLinear.cc \
	data/ElasticityImplicitLgDeformData3DQuadratic.cc \
	data/ElasticityImplicitLgDeformGravData2DLinear.cc \
	data/ElasticityImplicitLgDeformGravData2DQuadratic.cc \
	data/ElasticityImplicitLgDeformGravData3DLinear.cc \
	data/ElasticityImplicitLgDeformGravData3DQuadratic.cc \
	data/QuadratureData.cc \
	data/QuadratureData1Din2DLinear.cc \
	data/QuadratureData1Din2DQuadratic.cc \
	data/QuadratureData1Din3DLinear.cc \
	data/QuadratureData1Din3DQuadratic.cc \
	data/QuadratureData2DLinear.cc \
	data/QuadratureData2DQuadratic.cc \
	data/QuadratureData2Din3DLinearXYZ.cc \
	data/QuadratureData2Din3DLinearXY.cc \
	data/QuadratureData2Din3DLinearYZ.cc \
	data/QuadratureData2Din3DLinearXZ.cc \
	data/QuadratureData2Din3DQuadratic.cc \
	data/QuadratureData3DLinear.cc \
	data/QuadratureData3DQuadratic.cc

noinst_HEADERS += \
	data/CellGeomData.hh \
	data/GeomDataLine2D.hh \
	data/GeomDataLine3D.hh \
	data/GeomDataTri2D.hh \
	data/GeomDataTri3D.hh \
	data/GeomDataQuad2D.hh \
	data/GeomDataQuad3D.hh \
	data/GeomDataTet3D.hh \
	data/GeomDataHex3D.hh \
	data/IntegratorData.hh \
	data/ElasticityExplicitData.hh \
	data/ElasticityExplicitData2DLinear.hh \
	data/ElasticityExplicitData2DQuadratic.hh \
	data/ElasticityExplicitData3DLinear.hh \
	data/ElasticityExplicitData3DQuadratic.hh \
	data/ElasticityExplicitGravData2DLinear.hh \
	data/ElasticityExplicitGravData2DQuadratic.hh \
	data/ElasticityExplicitGravData3DLinear.hh \
	data/ElasticityExplicitGravData3DQuadratic.hh \
	data/ElasticityImplicitData2DLinear.hh \
	data/ElasticityImplicitData2DQuadratic.hh \
	data/ElasticityImplicitData3DLinear.hh \
	data/ElasticityImplicitData3DQuadratic.hh \
	data/ElasticityImplicitGravData2DLinear.hh \
	data/ElasticityImplicitGravData2DQuadratic.hh \
	data/ElasticityImplicitGravData3DLinear.hh \
	data/ElasticityImplicitGravData3DQuadratic.hh \
	data/ElasticityExplicitLgDeformData2DLinear.hh \
	data/ElasticityExplicitLgDeformData2DQuadratic.hh \
	data/ElasticityExplicitLgDeformData3DLinear.hh \
	data/ElasticityExplicitLgDeformData3DQuadratic.hh \
	data/ElasticityExplicitLgDeformGravData2DLinear.hh \
	data/ElasticityExplicitLgDeformGravData2DQuadratic.hh \
	data/ElasticityExplicitLgDeformGravData3DLinear.hh \
	data/ElasticityExplicitLgDeformGravData3DQuadratic.hh \
	data/ElasticityImplicitLgDeformData2DLinear.hh \
	data/ElasticityImplicitLgDeformData2DQuadratic.hh \
	data/ElasticityImplicitLgDeformData3DLinear.hh \
	data/ElasticityImplicitLgDeformData3DQuadratic.hh \
	data/ElasticityImplicitLgDeformGravData2DLinear.hh \
	data/ElasticityImplicitLgDeformGravData2DQuadratic.hh \
	data/ElasticityImplicitLgDeformGravData3DLinear.hh \
	data/ElasticityImplicitLgDeformGravData3DQuadratic.hh \
	data/QuadratureData.hh \
	data/QuadratureData1Din2DLinear.hh \
	data/QuadratureData1Din2DQuadratic.hh \
	data/QuadratureData1Din3DLinear.hh \
	data/QuadratureData1Din3DQuadratic.hh \
	data/QuadratureData2DLinear.hh \
	data/QuadratureData2DQuadratic.hh \
	data/QuadratureData2Din3DLinearXYZ.hh \
	data/QuadratureData2Din3DLinearXY.hh \
	data/QuadratureData2Din3DLinearYZ.hh \
	data/QuadratureData2Din3DLinearXZ.hh \
	data/QuadratureData2Din3DQuadratic.hh \
	data/QuadratureData3DLinear.hh \
	data/QuadratureData3DQuadratic.hh

AM_CPPFLAGS += \
	$(PYTHON_EGG_CPPFLAGS) -I$(PYTHON_INCDIR) \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES)

testfeassemble_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testfeassemble_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testfeassemble_LDADD += -lnetcdf_c++ -lnetcdf
endif


leakcheck: testfeassemble
	valgrind --log-file=valgrind_feassemble.log --leak-check=full --suppressions=$(top_srcdir)/share/valgrind-python.supp .libs/testfeassemble


# End of file 
