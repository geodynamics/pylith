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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmeshio.py
dist_check_SCRIPTS = testmeshio.py 

noinst_PYTHON = \
	TestMeshIOAscii.py \
	TestMeshIOCubit.py \
	TestMeshIOLagrit.py \
	TestVertexFilterVecNorm.py \
	TestCellFilterAvg.py \
	TestOutputManagerMesh.py \
	TestOutputManagerSubMesh.py \
	TestOutputSolnSubset.py \
	TestOutputSolnPoints.py \
	TestDataWriterVTK.py \
	TestDataWriterHDF5.py \
	TestDataWriterHDF5Ext.py \
	TestSingleOutput.py \
	TestXdmf.py


if ENABLE_HDF5
  TESTS += testhdf5.py
  dist_check_SCRIPTS += testhdf5.py
endif	

if ENABLE_CUBIT
  TESTS += testcubit.py
  dist_check_SCRIPTS += testcubit.py
endif	

noinst_tmp = \
	output_sub_t20.vtk

CLEANFILES = $(noinst_tmp)

# End of file 
