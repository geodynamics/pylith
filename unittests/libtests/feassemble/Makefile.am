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
# Copyright (c) 2010-2011 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

subpackage = feassemble
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfeassemble

check_PROGRAMS = testfeassemble

# Primary source files
testfeassemble_SOURCES = \
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
	TestGeometryQuad3D.cc \
	TestGeometryTet3D.cc \
	TestGeometryHex3D.cc \
	TestQuadratureRefCell.cc \
	TestQuadratureEngine.cc \
	TestQuadrature0D.cc \
	TestQuadrature1D.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature1Din3D.cc \
	TestQuadrature2D.cc \
	TestQuadrature2Din3D.cc \
	TestQuadrature3D.cc \
	TestQuadrature.cc \
	TestIntegrator.cc \
	TestIntegratorElasticity.cc \
	TestElasticityExplicit.cc \
	TestElasticityExplicit1DLinear.cc \
	TestElasticityExplicit1DQuadratic.cc \
	TestElasticityExplicit2DLinear.cc \
	TestElasticityExplicit2DQuadratic.cc \
	TestElasticityExplicit3DLinear.cc \
	TestElasticityExplicit3DQuadratic.cc \
	TestElasticityExplicitGrav1DLinear.cc \
	TestElasticityExplicitGrav1DQuadratic.cc \
	TestElasticityExplicitGrav2DLinear.cc \
	TestElasticityExplicitGrav2DQuadratic.cc \
	TestElasticityExplicitGrav3DLinear.cc \
	TestElasticityExplicitGrav3DQuadratic.cc \
	TestElasticityExplicitTri3.cc \
	TestElasticityExplicitTet4.cc \
	TestElasticityImplicit.cc \
	TestElasticityImplicit1DLinear.cc \
	TestElasticityImplicit1DQuadratic.cc \
	TestElasticityImplicit2DLinear.cc \
	TestElasticityImplicit2DQuadratic.cc \
	TestElasticityImplicit3DLinear.cc \
	TestElasticityImplicit3DQuadratic.cc \
	TestElasticityImplicitGrav1DLinear.cc \
	TestElasticityImplicitGrav1DQuadratic.cc \
	TestElasticityImplicitGrav2DLinear.cc \
	TestElasticityImplicitGrav2DQuadratic.cc \
	TestElasticityImplicitGrav3DLinear.cc \
	TestElasticityImplicitGrav3DQuadratic.cc \
	TestIntegratorElasticityLgDeform.cc \
	TestElasticityExplicitLgDeform.cc \
	TestElasticityExplicitLgDeform1DLinear.cc \
	TestElasticityExplicitLgDeform1DQuadratic.cc \
	TestElasticityExplicitLgDeform2DLinear.cc \
	TestElasticityExplicitLgDeform2DQuadratic.cc \
	TestElasticityExplicitLgDeform3DLinear.cc \
	TestElasticityExplicitLgDeform3DQuadratic.cc \
	TestElasticityExplicitLgDeformGrav1DLinear.cc \
	TestElasticityExplicitLgDeformGrav1DQuadratic.cc \
	TestElasticityExplicitLgDeformGrav2DLinear.cc \
	TestElasticityExplicitLgDeformGrav2DQuadratic.cc \
	TestElasticityExplicitLgDeformGrav3DLinear.cc \
	TestElasticityExplicitLgDeformGrav3DQuadratic.cc \
	TestElasticityImplicitLgDeform.cc \
	TestElasticityImplicitLgDeform1DLinear.cc \
	TestElasticityImplicitLgDeform1DQuadratic.cc \
	TestElasticityImplicitLgDeform2DLinear.cc \
	TestElasticityImplicitLgDeform2DQuadratic.cc \
	TestElasticityImplicitLgDeform3DLinear.cc \
	TestElasticityImplicitLgDeform3DQuadratic.cc \
	TestElasticityImplicitLgDeformGrav1DLinear.cc \
	TestElasticityImplicitLgDeformGrav1DQuadratic.cc \
	TestElasticityImplicitLgDeformGrav2DLinear.cc \
	TestElasticityImplicitLgDeformGrav2DQuadratic.cc \
	TestElasticityImplicitLgDeformGrav3DLinear.cc \
	TestElasticityImplicitLgDeformGrav3DQuadratic.cc \
	test_feassemble.cc



noinst_HEADERS = \
	TestCellGeometry.hh \
	TestElasticityExplicit.hh \
	TestElasticityExplicit1DLinear.hh \
	TestElasticityExplicit1DQuadratic.hh \
	TestElasticityExplicit2DLinear.hh \
	TestElasticityExplicit2DQuadratic.hh \
	TestElasticityExplicit3DLinear.hh \
	TestElasticityExplicit3DQuadratic.hh \
	TestElasticityExplicitGrav1DLinear.hh \
	TestElasticityExplicitGrav1DQuadratic.hh \
	TestElasticityExplicitGrav2DLinear.hh \
	TestElasticityExplicitGrav2DQuadratic.hh \
	TestElasticityExplicitGrav3DLinear.hh \
	TestElasticityExplicitGrav3DQuadratic.hh \
	TestElasticityExplicitTri3.hh \
	TestElasticityExplicitTet4.hh \
	TestElasticityImplicit.hh \
	TestElasticityImplicit1DLinear.hh \
	TestElasticityImplicit1DQuadratic.hh \
	TestElasticityImplicit2DLinear.hh \
	TestElasticityImplicit2DQuadratic.hh \
	TestElasticityImplicit3DLinear.hh \
	TestElasticityImplicit3DQuadratic.hh \
	TestElasticityImplicitGrav1DLinear.hh \
	TestElasticityImplicitGrav1DQuadratic.hh \
	TestElasticityImplicitGrav2DLinear.hh \
	TestElasticityImplicitGrav2DQuadratic.hh \
	TestElasticityImplicitGrav3DLinear.hh \
	TestElasticityImplicitGrav3DQuadratic.hh \
	TestIntegratorElasticityLgDeform.hh \
	TestElasticityExplicitLgDeform.hh \
	TestElasticityExplicitLgDeform1DLinear.hh \
	TestElasticityExplicitLgDeform1DQuadratic.hh \
	TestElasticityExplicitLgDeform2DLinear.hh \
	TestElasticityExplicitLgDeform2DQuadratic.hh \
	TestElasticityExplicitLgDeform3DLinear.hh \
	TestElasticityExplicitLgDeform3DQuadratic.hh \
	TestElasticityExplicitLgDeformGrav1DLinear.hh \
	TestElasticityExplicitLgDeformGrav1DQuadratic.hh \
	TestElasticityExplicitLgDeformGrav2DLinear.hh \
	TestElasticityExplicitLgDeformGrav2DQuadratic.hh \
	TestElasticityExplicitLgDeformGrav3DLinear.hh \
	TestElasticityExplicitLgDeformGrav3DQuadratic.hh \
	TestElasticityImplicitLgDeform.hh \
	TestElasticityImplicitLgDeform1DLinear.hh \
	TestElasticityImplicitLgDeform1DQuadratic.hh \
	TestElasticityImplicitLgDeform2DLinear.hh \
	TestElasticityImplicitLgDeform2DQuadratic.hh \
	TestElasticityImplicitLgDeform3DLinear.hh \
	TestElasticityImplicitLgDeform3DQuadratic.hh \
	TestElasticityImplicitLgDeformGrav1DLinear.hh \
	TestElasticityImplicitLgDeformGrav1DQuadratic.hh \
	TestElasticityImplicitLgDeformGrav2DLinear.hh \
	TestElasticityImplicitLgDeformGrav2DQuadratic.hh \
	TestElasticityImplicitLgDeformGrav3DLinear.hh \
	TestElasticityImplicitLgDeformGrav3DQuadratic.hh \
	TestGeometryPoint1D.hh \
	TestGeometryPoint2D.hh \
	TestGeometryPoint3D.hh \
	TestGeometryLine1D.hh \
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
	TestQuadrature0D.hh \
	TestQuadrature1D.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature1Din3D.hh \
	TestQuadrature2D.hh \
	TestQuadrature2Din3D.hh \
	TestQuadrature3D.hh

# Source files associated with testing data
testfeassemble_SOURCES += \
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
	data/GeomDataQuad3D.cc \
	data/GeomDataTet3D.cc \
	data/GeomDataHex3D.cc \
	data/IntegratorData.cc \
	data/ElasticityExplicitData.cc \
	data/ElasticityExplicitData1DLinear.cc \
	data/ElasticityExplicitData1DQuadratic.cc \
	data/ElasticityExplicitData2DLinear.cc \
	data/ElasticityExplicitData2DQuadratic.cc \
	data/ElasticityExplicitData3DLinear.cc \
	data/ElasticityExplicitData3DQuadratic.cc \
	data/ElasticityExplicitGravData1DLinear.cc \
	data/ElasticityExplicitGravData1DQuadratic.cc \
	data/ElasticityExplicitGravData2DLinear.cc \
	data/ElasticityExplicitGravData2DQuadratic.cc \
	data/ElasticityExplicitGravData3DLinear.cc \
	data/ElasticityExplicitGravData3DQuadratic.cc \
	data/ElasticityImplicitData1DLinear.cc \
	data/ElasticityImplicitData1DQuadratic.cc \
	data/ElasticityImplicitData2DLinear.cc \
	data/ElasticityImplicitData2DQuadratic.cc \
	data/ElasticityImplicitData3DLinear.cc \
	data/ElasticityImplicitData3DQuadratic.cc \
	data/ElasticityImplicitGravData1DLinear.cc \
	data/ElasticityImplicitGravData1DQuadratic.cc \
	data/ElasticityImplicitGravData2DLinear.cc \
	data/ElasticityImplicitGravData2DQuadratic.cc \
	data/ElasticityImplicitGravData3DLinear.cc \
	data/ElasticityImplicitGravData3DQuadratic.cc \
	data/ElasticityExplicitLgDeformData1DLinear.cc \
	data/ElasticityExplicitLgDeformData1DQuadratic.cc \
	data/ElasticityExplicitLgDeformData2DLinear.cc \
	data/ElasticityExplicitLgDeformData2DQuadratic.cc \
	data/ElasticityExplicitLgDeformData3DLinear.cc \
	data/ElasticityExplicitLgDeformData3DQuadratic.cc \
	data/ElasticityExplicitLgDeformGravData1DLinear.cc \
	data/ElasticityExplicitLgDeformGravData1DQuadratic.cc \
	data/ElasticityExplicitLgDeformGravData2DLinear.cc \
	data/ElasticityExplicitLgDeformGravData2DQuadratic.cc \
	data/ElasticityExplicitLgDeformGravData3DLinear.cc \
	data/ElasticityExplicitLgDeformGravData3DQuadratic.cc \
	data/ElasticityImplicitLgDeformData1DLinear.cc \
	data/ElasticityImplicitLgDeformData1DQuadratic.cc \
	data/ElasticityImplicitLgDeformData2DLinear.cc \
	data/ElasticityImplicitLgDeformData2DQuadratic.cc \
	data/ElasticityImplicitLgDeformData3DLinear.cc \
	data/ElasticityImplicitLgDeformData3DQuadratic.cc \
	data/ElasticityImplicitLgDeformGravData1DLinear.cc \
	data/ElasticityImplicitLgDeformGravData1DQuadratic.cc \
	data/ElasticityImplicitLgDeformGravData2DLinear.cc \
	data/ElasticityImplicitLgDeformGravData2DQuadratic.cc \
	data/ElasticityImplicitLgDeformGravData3DLinear.cc \
	data/ElasticityImplicitLgDeformGravData3DQuadratic.cc \
	data/QuadratureData.cc \
	data/QuadratureData1DLinear.cc \
	data/QuadratureData1DQuadratic.cc \
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
	data/GeomDataPoint1D.hh \
	data/GeomDataPoint2D.hh \
	data/GeomDataPoint3D.hh \
	data/GeomDataLine1D.hh \
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
	data/ElasticityExplicitData1DLinear.hh \
	data/ElasticityExplicitData1DQuadratic.hh \
	data/ElasticityExplicitData2DLinear.hh \
	data/ElasticityExplicitData2DQuadratic.hh \
	data/ElasticityExplicitData3DLinear.hh \
	data/ElasticityExplicitData3DQuadratic.hh \
	data/ElasticityExplicitGravData1DLinear.hh \
	data/ElasticityExplicitGravData1DQuadratic.hh \
	data/ElasticityExplicitGravData2DLinear.hh \
	data/ElasticityExplicitGravData2DQuadratic.hh \
	data/ElasticityExplicitGravData3DLinear.hh \
	data/ElasticityExplicitGravData3DQuadratic.hh \
	data/ElasticityImplicitData1DLinear.hh \
	data/ElasticityImplicitData1DQuadratic.hh \
	data/ElasticityImplicitData2DLinear.hh \
	data/ElasticityImplicitData2DQuadratic.hh \
	data/ElasticityImplicitData3DLinear.hh \
	data/ElasticityImplicitData3DQuadratic.hh \
	data/ElasticityImplicitGravData1DLinear.hh \
	data/ElasticityImplicitGravData1DQuadratic.hh \
	data/ElasticityImplicitGravData2DLinear.hh \
	data/ElasticityImplicitGravData2DQuadratic.hh \
	data/ElasticityImplicitGravData3DLinear.hh \
	data/ElasticityImplicitGravData3DQuadratic.hh \
	data/ElasticityExplicitLgDeformData1DLinear.hh \
	data/ElasticityExplicitLgDeformData1DQuadratic.hh \
	data/ElasticityExplicitLgDeformData2DLinear.hh \
	data/ElasticityExplicitLgDeformData2DQuadratic.hh \
	data/ElasticityExplicitLgDeformData3DLinear.hh \
	data/ElasticityExplicitLgDeformData3DQuadratic.hh \
	data/ElasticityExplicitLgDeformGravData1DLinear.hh \
	data/ElasticityExplicitLgDeformGravData1DQuadratic.hh \
	data/ElasticityExplicitLgDeformGravData2DLinear.hh \
	data/ElasticityExplicitLgDeformGravData2DQuadratic.hh \
	data/ElasticityExplicitLgDeformGravData3DLinear.hh \
	data/ElasticityExplicitLgDeformGravData3DQuadratic.hh \
	data/ElasticityImplicitLgDeformData1DLinear.hh \
	data/ElasticityImplicitLgDeformData1DQuadratic.hh \
	data/ElasticityImplicitLgDeformData2DLinear.hh \
	data/ElasticityImplicitLgDeformData2DQuadratic.hh \
	data/ElasticityImplicitLgDeformData3DLinear.hh \
	data/ElasticityImplicitLgDeformData3DQuadratic.hh \
	data/ElasticityImplicitLgDeformGravData1DLinear.hh \
	data/ElasticityImplicitLgDeformGravData1DQuadratic.hh \
	data/ElasticityImplicitLgDeformGravData2DLinear.hh \
	data/ElasticityImplicitLgDeformGravData2DQuadratic.hh \
	data/ElasticityImplicitLgDeformGravData3DLinear.hh \
	data/ElasticityImplicitLgDeformGravData3DQuadratic.hh \
	data/QuadratureData.hh \
	data/QuadratureData1DLinear.hh \
	data/QuadratureData1DQuadratic.hh \
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

AM_CPPFLAGS = \
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


# End of file 
