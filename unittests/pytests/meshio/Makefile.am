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
dist_check_SCRIPTS = testdriver.py 

TESTS_ENVIRONMENT = $(PYTHON)

noinst_PYTHON = \
	TestMeshIOAscii.py \
	TestMeshIOCubit.py \
	TestMeshIOLagrit.py \
	TestVertexFilterVecNorm.py \
	TestCellFilterAvg.py \
	TestOutputManager.py \
	TestOutputSolnSubset.py

if ENABLE_CUBIT
#  TESTS += testdriver_cubit.py
#  dist_check_SCRIPTS += testdriver_cubit.py
endif	

noinst_tmp = \
	output_info.vtk \
	output_t20.vtk \
	output_sub_t20.vtk

CLEANFILES = $(noinst_tmp)

# End of file 
