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

subpackage = topology
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testtopology.py
dist_check_SCRIPTS = testtopology.py

TESTS_ENVIRONMENT = $(PYTHON)

#if ENABLE_TETGEN
#  TESTS += testgensimple.py
#  dist_check_SCRIPTS += testgensimple.py
#endif


noinst_PYTHON = \
	TestMesh.py \
	TestSubMesh.py \
	TestFieldBase.py \
	TestMeshField.py \
	TestMeshFields.py \
	TestSolutionFields.py \
	TestJacobian.py \
	TestMeshGenerator.py \
	TestMeshGenSimple.py \
	TestMeshImporter.py \
	TestRefineUniform.py


noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

CLEANFILES = $(noinst_tmp)


# End of file 
