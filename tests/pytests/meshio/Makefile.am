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
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmeshio.py testhdf5.py
dist_check_SCRIPTS = $(TESTS)

noinst_PYTHON = \
	TestMeshIOAscii.py \
	TestMeshIOCubit.py \
	TestMeshIOLagrit.py \
	TestVertexFilterVecNorm.py \
	TestCellFilterAvg.py \
	TestOutputManagerMesh.py \
	TestOutputManagerSubmesh.py \
	TestOutputSolnSubset.py \
	TestOutputSolnPoints.py \
	TestDataWriterVTK.py \
	TestDataWriterHDF5.py \
	TestDataWriterHDF5Ext.py \
	TestSingleOutput.py \
	TestXdmf.py


if ENABLE_CUBIT
  TESTS += testcubit.py
  dist_check_SCRIPTS += testcubit.py
endif

noinst_tmp = \
	output_sub_t20.vtk

CLEANFILES = $(noinst_tmp)

clean-local: clean-local-tmp
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.vtk *.dat *.dat.info *.h5 *.xmf

# End of file
