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

TESTS = testdriver.py

check_SCRIPTS = testdriver.py

noinst_PYTHON = \
	TestMesh.py \
	TestMeshGenerator.py \
	TestMeshImporter.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
