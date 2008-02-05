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

subpackage = meshio
include $(top_srcdir)/subpackage.am

SUBDIRS = data

TESTS = testdriver.py
check_SCRIPTS = testdriver.py 

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestMeshIOAscii.py \
	TestMeshIOCubit.py \
	TestMeshIOLagrit.py

if ENABLE_CUBIT
  TESTS += testdriver_cubit.py
  check_SCRIPTS += testdriver_cubit.py
endif	


# End of file 
