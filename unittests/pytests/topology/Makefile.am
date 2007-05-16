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

if ENABLE_TETGEN
  TESTS += testgensimple.py
  check_SCRIPTS += testgensimple.py
endif


noinst_PYTHON = \
	TestMesh.py \
	TestMeshGenerator.py \
	TestMeshGenSimple.py \
	TestMeshImporter.py

TESTS_ENVIRONMENT = $(PYTHON)


# End of file 
