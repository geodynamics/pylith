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

TESTS = testfeassemble

check_PROGRAMS = testfeassemble

# Primary source files
testfeassemble_SOURCES = \
	TestElasticityExplicit.cc \
	TestIntegrator.cc \
	TestQuadrature.cc \
	TestQuadrature1D.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature2D.cc \
	TestQuadrature2Din3D.cc \
	TestQuadrature3D.cc \
	test_feassemble.cc

noinst_HEADERS = \
	TestElasticityExplicit.hh \
	TestIntegrator.hh \
	TestQuadrature.hh \
	TestQuadrature1D.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature1Din3D.hh \
	TestQuadrature2D.hh \
	TestQuadrature2Din3D.hh

# Source files associated with testing data
testfeassemble_SOURCES += \
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
	$(top_builddir)/libsrc/libpylith.la

if ENABLE_CUBIT
  testfeassemble_LDADD += -lnetcdf_c++ -lnetcdf
endif


# End of file 
