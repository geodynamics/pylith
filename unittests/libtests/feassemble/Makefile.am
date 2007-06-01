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

subpackage = feassemble
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testfeassemble

check_PROGRAMS = testfeassemble

# Primary source files
testfeassemble_SOURCES = \
	TestElasticity.cc \
	TestElasticityExplicit.cc \
	TestElasticityExplicit1DLinear.cc \
	TestElasticityExplicit1DQuadratic.cc \
	TestElasticityExplicit2DLinear.cc \
	TestElasticityExplicit2DQuadratic.cc \
	TestElasticityExplicit3DLinear.cc \
	TestElasticityExplicit3DQuadratic.cc \
	TestElasticityImplicit.cc \
	TestElasticityImplicit1DLinear.cc \
	TestElasticityImplicit1DQuadratic.cc \
	TestElasticityImplicit2DLinear.cc \
	TestElasticityImplicit2DQuadratic.cc \
	TestElasticityImplicit3DLinear.cc \
	TestElasticityImplicit3DQuadratic.cc \
	TestIntegrator.cc \
	TestQuadrature.cc \
	TestQuadrature1D.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature2D.cc \
	TestQuadrature2Din3D.cc \
	TestQuadrature3D.cc \
	test_feassemble.cc

noinst_HEADERS = \
	TestElasticity.hh \
	TestElasticityExplicit.hh \
	TestElasticityExplicit1DLinear.hh \
	TestElasticityExplicit1DQuadratic.hh \
	TestElasticityExplicit2DLinear.hh \
	TestElasticityExplicit2DQuadratic.hh \
	TestElasticityExplicit3DLinear.hh \
	TestElasticityExplicit3DQuadratic.hh \
	TestElasticityImplicit.hh \
	TestElasticityImplicit1DLinear.hh \
	TestElasticityImplicit1DQuadratic.hh \
	TestElasticityImplicit2DLinear.hh \
	TestElasticityImplicit2DQuadratic.hh \
	TestElasticityImplicit3DLinear.hh \
	TestElasticityImplicit3DQuadratic.hh \
	TestIntegrator.hh \
	TestQuadrature.hh \
	TestQuadrature1D.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature1Din3D.hh \
	TestQuadrature2D.hh \
	TestQuadrature2Din3D.hh

# Source files associated with testing data
testfeassemble_SOURCES += \
	data/IntegratorData.cc \
	data/ElasticityExplicitData1DLinear.cc \
	data/ElasticityExplicitData1DQuadratic.cc \
	data/ElasticityExplicitData2DLinear.cc \
	data/ElasticityExplicitData2DQuadratic.cc \
	data/ElasticityExplicitData3DLinear.cc \
	data/ElasticityExplicitData3DQuadratic.cc \
	data/ElasticityImplicitData1DLinear.cc \
	data/ElasticityImplicitData1DQuadratic.cc \
	data/ElasticityImplicitData2DLinear.cc \
	data/ElasticityImplicitData2DQuadratic.cc \
	data/ElasticityImplicitData3DLinear.cc \
	data/ElasticityImplicitData3DQuadratic.cc \
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
	data/IntegratorData.hh \
	data/ElasticityExplicitData1DLinear.hh \
	data/ElasticityExplicitData1DQuadratic.hh \
	data/ElasticityExplicitData2DLinear.hh \
	data/ElasticityExplicitData2DQuadratic.hh \
	data/ElasticityExplicitData3DLinear.hh \
	data/ElasticityExplicitData3DQuadratic.hh \
	data/ElasticityImplicitData1DLinear.hh \
	data/ElasticityImplicitData1DQuadratic.hh \
	data/ElasticityImplicitData2DLinear.hh \
	data/ElasticityImplicitData2DQuadratic.hh \
	data/ElasticityImplicitData3DLinear.hh \
	data/ElasticityImplicitData3DQuadratic.hh \
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

testfeassemble_LDFLAGS = $(PETSC_LIB) $(PYTHON_BLDLIBRARY)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfeassemble_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la \
	-lspatialdata

if ENABLE_CUBIT
  testfeassemble_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
