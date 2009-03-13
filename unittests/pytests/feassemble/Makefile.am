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

TESTS = testfeassemble.py
dist_check_SCRIPTS = testfeassemble.py

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestFIATLagrange.py \
	TestFIATSimplex.py \
	TestCellGeometry.py \
	TestMeshQuadrature.py \
	TestElasticityExplicit.py \
	TestElasticityImplicit.py


# End of file 
