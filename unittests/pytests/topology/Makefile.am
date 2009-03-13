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

subpackage = topology
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testdriver.py
dist_check_SCRIPTS = testdriver.py

TESTS_ENVIRONMENT = $(PYTHON)

#if ENABLE_TETGEN
#  TESTS += testgensimple.py
#  dist_check_SCRIPTS += testgensimple.py
#endif


noinst_PYTHON = \
	TestMesh.py \
	TestSubMesh.py \
	TestFieldsManager.py \
	TestMeshGenerator.py \
	TestMeshGenSimple.py \
	TestMeshImporter.py


noinst_tmp = \
	jacobian.mat \
	jacobian.mat.info

CLEANFILES = $(noinst_tmp)


# End of file 
