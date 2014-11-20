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

TESTS = testpylith.py

dist_check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestTet4.py \
	TestAxialDisp.py \
	axialdisp_gendb.py \
	axialdisp_soln.py \
	TestShearDisp.py \
	TestShearDispNoSlip.py \
	TestShearDispNoSlipRefine.py \
	sheardisp_gendb.py \
	sheardisp_soln.py \
	TestSlipOneFault.py \
	sliponefault_soln.py \
	TestSlipTwoFaults.py \
	sliptwofaults_soln.py \
	TestFaultsIntersect.py


dist_noinst_DATA = \
	geometry.jou \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	axialdisp.cfg \
	sheardisp.cfg \
	sheardispnoslip.cfg \
	sheardispnosliprefine.cfg \
	sliponefault.cfg \
	points.txt \
	sliptwofaults.cfg \
	faultsintersect.cfg

noinst_TMP = \
	axial_dispx.spatialdb \
	axial_dispy.spatialdb \
	axial_dispz.spatialdb \
	shear_dispx.spatialdb \
	shear_dispy.spatialdb \
	shear_dispz.spatialdb


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/3d/tet4
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.h5 *.xmf *.pyc *.dat *.dat.info


# End of file 
