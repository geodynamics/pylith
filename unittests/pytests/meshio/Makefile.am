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

noinst_PYTHON = \
	TestMeshIOAscii.py \
	TestMeshIOCubit.py \
	TestMeshIOLagrit.py

if ENABLE_CUBIT
  TESTS += testdriver_cubit.py
  check_SCRIPTS += testdriver_cubit.py
endif	

TESTS_ENVIRONMENT = $(PYTHON)




# End of file 
