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
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

dist_noinst_DATA = \
	README.md \
	pylithapp.cfg \
	greensfns.cfg \
	mat_elastic.cfg \
	mat_viscoelastic.cfg \
	solver_algebraicmultigrid.cfg \
	solver_fieldsplit.cfg \
	step01.cfg \
	step02.cfg \
	step03.cfg \
	step04.cfg \
	step05.cfg \
	step06.cfg \
	step07a.cfg \
	step07b.cfg \
	step08a.cfg \
	step08b.cfg \
	step08c.cfg \
	cgps_sites.txt \
	make_synthetic_gpsdisp.cfg \
	slip_invert.cfg \
	make_synthetic_gpsdisp.py \
	slip_invert.py \
	mesh/bc.jou \
	mesh/cas_contours_dep.in.txt.gz \
	mesh/coordsys.py \
	mesh/generate_surfjou.py \
	mesh/geometry.jou \
	mesh/__init__.py \
	mesh/mesh_hex.jou \
	mesh/mesh_tet.jou \
	spatialdb/mat_crust_elastic.spatialdb \
	spatialdb/mat_mantle_elastic.spatialdb \
	spatialdb/mat_slab_elastic.spatialdb \
	spatialdb/mat_wedge_elastic.spatialdb \
	spatialdb/mat_viscosity.spatialdb \
	spatialdb/mat_initial_stress_grav.spatialdb \
	spatialdb/fault_slabtop_coseismic.spatialdb \
	spatialdb/fault_slabtop_creep.spatialdb \
	spatialdb/fault_slabtop_slipweakening.spatialdb \
	spatialdb/fault_slabtop_tractions.spatialdb \
	spatialdb/generate_initial_stress.py \
	spatialdb/generate_slowslip.cfg \
	spatialdb/generate_slowslip.py \
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
	viz/cascadia_tm.vtk


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/examples/3d/subduction
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
