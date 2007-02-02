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
	TestIntegrator.cc \
	TestIntegratorInertia.cc \
	TestIntegratorInertia1D.cc \
	TestIntegratorInertia2D.cc \
	TestIntegratorInertia3D.cc \
	TestQuadrature.cc \
	TestQuadrature1D.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature2D.cc \
	TestQuadrature2Din3D.cc \
	TestQuadrature3D.cc \
	test_feassemble.cc

noinst_HEADERS = \
	TestIntegrator.hh \
	TestIntegratorInertia.hh \
	TestIntegratorInertia1D.hh \
	TestIntegratorInertia2D.hh \
	TestIntegratorInertia3D.hh \
	TestQuadrature.hh \
	TestQuadrature1D.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature1Din3D.hh \
	TestQuadrature2D.hh \
	TestQuadrature2Din3D.hh

# Source files associated with testing data
testfeassemble_SOURCES += \
	data/IntegratorData.cc \
	data/IntegratorDataInertia1DLinear.cc \
	data/IntegratorDataInertia1DQuadratic.cc \
	data/IntegratorDataInertia2Din3DOne.cc \
	data/IntegratorDataInertia2Din3DTwo.cc \
	data/IntegratorDataInertia2Din3DThree.cc \
	data/IntegratorDataInertia3DLinear.cc \
	data/IntegratorDataInertia3DQuadratic.cc \
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
	data/IntegratorDataInertia1DLinear.hh \
	data/IntegratorDataInertia1DQuadratic.hh \
	data/IntegratorDataInertia2Din3DOne.hh \
	data/IntegratorDataInertia2Din3DTwo.hh \
	data/IntegratorDataInertia2Din3DThree.hh \
	data/IntegratorDataInertia3DLinear.hh \
	data/IntegratorDataInertia3DQuadratic.hh \
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

testfeassemble_LDFLAGS = $(PETSC_LIB)

AM_CPPFLAGS = $(PETSC_SIEVE_FLAGS) $(PETSC_INCLUDE)

testfeassemble_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/libpylith.la

# End of file 
