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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testtopology.py
dist_check_SCRIPTS = testtopology.py

noinst_PYTHON = \
	TestMesh.py \
	TestSubMesh.py \
	TestFieldBase.py \
	TestMeshField.py \
	TestMeshFields.py \
	TestSolutionFields.py \
	TestJacobian.py \
	TestMeshGenerator.py \
	TestMeshImporter.py \
	TestRefineUniform.py


noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

CLEANFILES = $(noinst_tmp)


# End of file 
