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

TESTS = testmeshio.py
dist_check_SCRIPTS = testmeshio.py 

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
  TESTS += testcubit.py
  dist_check_SCRIPTS += testcubit.py
endif	

noinst_tmp = \
	output_info.vtk \
	output_t20.vtk \
	output_sub_t20.vtk

CLEANFILES = $(noinst_tmp)

# End of file 
