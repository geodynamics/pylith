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
	test_feassemble.cc \
	data/QuadratureData.cc \
	data/QuadratureData1DLinear.cc \
	data/QuadratureData1DQuadratic.cc

noinst_HEADERS = \
	TestQuadrature.hh \
	TestQuadrature1D.hh
	data/QuadratureData.hh \
	data/QuadratureData1DLinear.hh \
	data/QuadratureData1DQuadratic.hh

testfeassemble_LDFLAGS = $(PETSC_LIB)

INCLUDES += $(PETSC_INCLUDE)

testfeassemble_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/feassemble/libpylithfeassemble.la

# version
# $Id$

# End of file 
