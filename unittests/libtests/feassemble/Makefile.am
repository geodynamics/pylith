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

testfeassemble_SOURCES = \
	TestQuadrature.cc \
	TestQuadrature1D.cc \
	TestQuadrature1Din2D.cc \
	TestQuadrature2D.cc \
	test_feassemble.cc \
	data/QuadratureData.cc \
	data/QuadratureData1DLinear.cc \
	data/QuadratureData1DQuadratic.cc \
	data/QuadratureData1Din2DLinear.cc \
	data/QuadratureData1Din2DQuadratic.cc \
	data/QuadratureData1Din3DLinear.cc \
	data/QuadratureData1Din3DQuadratic.cc \
	data/QuadratureData2DLinear.cc \
	data/QuadratureData2DQuadratic.cc

noinst_HEADERS = \
	TestQuadrature.hh \
	TestQuadrature1D.hh \
	TestQuadrature1Din2D.hh \
	TestQuadrature2D.hh \
	data/QuadratureData.hh \
	data/QuadratureData1DLinear.hh \
	data/QuadratureData1DQuadratic.hh \
	data/QuadratureData1Din2DLinear.hh \
	data/QuadratureData1Din2DQuadratic.hh \
	data/QuadratureData1Din3DLinear.hh \
	data/QuadratureData1Din3DQuadratic.hh \
	data/QuadratureData2DLinear.hh \
	data/QuadratureData2DQuadratic.hh 

testfeassemble_LDFLAGS = $(PETSC_LIB)

INCLUDES += $(PETSC_INCLUDE)

testfeassemble_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/feassemble/libpylithfeassemble.la

# version
# $Id$

# End of file 
