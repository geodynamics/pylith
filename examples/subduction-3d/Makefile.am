# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

dist_noinst_DATA = \
	README.md \
	cubit_bc.jou \
	cubit_geometry.jou \
	cubit_hex.jou \
	cubit_tet.jou \
	fault_slabtop_coseismic.spatialdb \
	fault_slabtop_creep.spatialdb \
	fault_slabtop_slipweakening.spatialdb \
	fault_slabtop_tractions.spatialdb \
	mat_elastic.cfg \
	mat_crust_elastic.spatialdb \
	mat_mantle_elastic.spatialdb \
	mat_slab_elastic.spatialdb \
	mat_wedge_elastic.spatialdb \
	mat_viscoelastic.cfg \
	mat_viscosity.spatialdb \
	mat_elastic_incompressible.cfg \
	mat_initial_stress_grav.spatialdb \
	pylithapp.cfg \
	step01_axialdisp.cfg \
	step02_coseismic.cfg \
	step03_interseismic.cfg \
	step04_eqcycle.cfg \
	step06_slowslip.cfg \
	step07a_leftlateral.cfg \
	step07b_reverse.cfg \
	step08a_gravity_refstate.cfg \
	step08b_gravity_incompressible.cfg \
	step08c_gravity_viscoelastic.cfg \
	cas_contours_dep.in.txt.gz \
	utils/coordsys.py \
	utils/generate_initial_stress.py \
	utils/generate_slowslip.py \
	utils/generate_surfjou.py \
	utils/make_synthetic_gpsdisp.py \
	utils/slip_invert.py \
	viz/plot_dispvec.py \
	viz/plot_dispwarp.py \
	viz/plot_faultdir.py \
	viz/plot_faultslip.py \
	viz/plot_faultslip-reltract.py \
	viz/plot_matprops.py \
	viz/plot_mesh.py \
	viz/plot_inversion_misfit.py \
	viz/get_coast.sh \
	viz/seg2vtk.py \
	viz/cascadia_tm.vtk \
	input/README.md \
	scratch/README.md


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/subduction-3d
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.pyc output/*


# End of file 
