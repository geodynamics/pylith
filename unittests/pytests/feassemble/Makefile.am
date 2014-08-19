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

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testfeassemble.py
dist_check_SCRIPTS = testfeassemble.py

noinst_PYTHON = \
	TestFIATLagrange.py \
	TestFIATSimplex.py \
	TestCellGeometry.py \
	TestMeshQuadrature.py \
	TestElasticityExplicit.py \
	TestElasticityImplicit.py \
	TestElasticityExplicitLgDeform.py \
	TestElasticityImplicitLgDeform.py


# End of file 
