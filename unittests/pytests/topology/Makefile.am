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
check_SCRIPTS = testdriver.py

TESTS_ENVIRONMENT = $(PYTHON)

if ENABLE_TETGEN
  TESTS += testgensimple.py
  check_SCRIPTS += testgensimple.py
endif


noinst_PYTHON = \
	TestMesh.py \
	TestMeshGenerator.py \
	TestMeshGenSimple.py \
	TestMeshImporter.py


# End of file 
