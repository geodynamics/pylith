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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

TESTS = testpylith.py

check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestQuad4.py \
	TestAxialDisp.py \
	axialdisp_soln.py \
	axialdisp_gendb.py \
	TestShearDisp.py \
	sheardisp_soln.py \
	sheardisp_gendb.py \
	TestDislocation.py \
	TestDislocation2.py \
	dislocation_soln.py \
	dislocation2_soln.py \
	TestLgDeformRigidBody.py \
	rigidbody_soln.py \
	rigidbody_gendb.py \
	TestLgDeformTraction.py \
	lgdeformtraction_soln.py \
	TestFrictionCompression.py \
	TestFrictionOpening.py \
	TestFrictionShearSliding.py \
	TestFrictionShearStick.py \
	TestSlipWeakeningCompression.py \
	TestSlipWeakeningOpening.py \
	TestSlipWeakeningShearSliding.py \
	TestSlipWeakeningShearStick.py \
	friction_compression_soln.py \
	friction_shear_stick_soln.py \
	friction_shear_sliding_soln.py \
	friction_opening_soln.py \
	slipweakening_compression_soln.py \
	slipweakening_shear_stick_soln.py \
	slipweakening_shear_sliding_soln.py \
	slipweakening_opening_soln.py \
	testpylith.py

dist_noinst_DATA = \
	geometry.jou \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	axialdisp.cfg \
	sheardisp.cfg \
	dislocation.cfg \
	dislocation2.cfg \
	lgdeformrigidbody.cfg \
	lgdeformtraction.cfg \
	friction_compression.cfg \
	friction_shear_stick.cfg \
	friction_shear_sliding.cfg \
	friction_opening.cfg \
	slipweakening_compression.cfg \
	slipweakening_shear_stick.cfg \
	slipweakening_shear_sliding.cfg \
	slipweakening_opening.cfg

noinst_TMP = \
	axial_disp.spatialdb \
	axialdisp_t0000000.vtk \
	axialdisp-elastic_info.vtk \
	axialdisp-elastic_t0000000.vtk \
	shear_disp.spatialdb \
	sheardisp_t0000000.vtk \
	sheardisp-elastic_info.vtk \
	sheardisp-elastic_t0000000.vtk \
	dislocation_t0000000.vtk \
	dislocation-elastic_info.vtk \
	dislocation-elastic_t0000000.vtk \
	dislocation-fault_info.vtk \
	dislocation-fault_t0000000.vtk \
	rigidbody_disp.spatialdb \
	lgdeformrigidbody_t0000000.vtk \
	lgdeformrigidbody-elastic_info.vtk \
	lgdeformrigidbody-elastic_t0000000.vtk \
	lgdeformtraction_t0000000.vtk \
	lgdeformtraction-elastic_info.vtk \
	lgdeformtraction-elastic_t0000000.vtk \
	friction_compression_t0000000.vtk \
	friction_compression-elastic_info.vtk \
	friction_compression-elastic_t0000000.vtk \
	friction_shear_stick_t0000000.vtk \
	friction_shear_stick-elastic_info.vtk \
	friction_shear_stick-elastic_t0000000.vtk \
	friction_shear_sliding_t0000000.vtk \
	friction_shear_sliding-elastic_info.vtk \
	friction_shear_sliding-elastic_t0000000.vtk \
	friction_opening_t0000000.vtk \
	friction_opening-elastic_info.vtk \
	friction_opening-elastic_t0000000.vtk \
	slipweakening_compression_t0000000.vtk \
	slipweakening_compression-elastic_info.vtk \
	slipweakening_compression-elastic_t0000000.vtk \
	slipweakening_shear_stick_t0000000.vtk \
	slipweakening_shear_stick-elastic_info.vtk \
	slipweakening_shear_stick-elastic_t0000000.vtk \
	slipweakening_shear_sliding_t0000000.vtk \
	slipweakening_shear_sliding-elastic_info.vtk \
	slipweakening_shear_sliding-elastic_t0000000.vtk \
	slipweakening_opening_t0000000.vtk \
	slipweakening_opening-elastic_info.vtk \
	slipweakening_opening-elastic_t0000000.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/2d/quad4
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc

# End of file 
