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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

TESTS = testpylith.py

check_SCRIPTS = testpylith.py

noinst_PYTHON = \
	TestAxialElasticIsotropic.py

noinst_DATA = \
	box.exo \
	matprops.spatialdb \
	matprops-maxwell.spatialdb \
	axial-small-tractions.cfg \
	axialelasticisotropic.cfg \
	axialgenmax1-tract-8cells.cfg \
	axialgenmax1a-tract-8cells.cfg \
	axialgenmax2-tract-8cells.cfg \
	axialgenmax3-tract-8cells.cfg \
	axialmaxwellisotropic.cfg \
	axialtractions-small.spatialdb \
	axialxyz.spatialdb \
	axialx.spatialdb \
	axialy.spatialdb \
	axialz.spatialdb \
	genmaxwell1.spatialdb \
	genmaxwell1a.spatialdb \
	genmaxwell2.spatialdb \
	genmaxwell3.spatialdb \
	pylithapp.cfg \
	shearelasticisotropic.cfg \
	shearplanestrain.cfg \
	shearmaxwellisotropic.cfg \
	shearxyz.spatialdb \
	twocells.mesh \
	twocellsb.mesh \
	threecells.mesh \
	fourcells.mesh \
	eightcells.mesh \
	shearelastic-small.cfg \
	shearelastic-twocells.cfg \
	shearelastic-twocellsb.cfg \
	shearelastic-threecells.cfg \
	shearelastic-fourcells.cfg \
	shearelastic-eightcells.cfg \
	shearxyz-small.spatialdb

noinst_TMP = \
	axialelasticisotropic_t0.vtk \
	axialelasticisotropic_t1.vtk \
	axialgenmax1-tract-8cells_t0000000.vtk \
	axialgenmax1-tract-8cells_t11045160000000.vtk \
	axialgenmax1-tract-8cells_t12623040000000.vtk \
	axialgenmax1-tract-8cells_t14200920000000.vtk \
	axialgenmax1-tract-8cells_t1577880000000.vtk \
	axialgenmax1-tract-8cells_t15778800000000.vtk \
	axialgenmax1-tract-8cells_t3155760000000.vtk \
	axialgenmax1-tract-8cells_t4733640000000.vtk \
	axialgenmax1-tract-8cells_t6311520000000.vtk \
	axialgenmax1-tract-8cells_t7889400000000.vtk \
	axialgenmax1-tract-8cells_t9467280000000.vtk \
	axialgenmax1a-tract-8cells_t0000000.vtk \
	axialgenmax1a-tract-8cells_t11045160000000.vtk \
	axialgenmax1a-tract-8cells_t12623040000000.vtk \
	axialgenmax1a-tract-8cells_t14200920000000.vtk \
	axialgenmax1a-tract-8cells_t1577880000000.vtk \
	axialgenmax1a-tract-8cells_t15778800000000.vtk \
	axialgenmax1a-tract-8cells_t3155760000000.vtk \
	axialgenmax1a-tract-8cells_t4733640000000.vtk \
	axialgenmax1a-tract-8cells_t6311520000000.vtk \
	axialgenmax1a-tract-8cells_t7889400000000.vtk \
	axialgenmax1a-tract-8cells_t9467280000000.vtk \
	axialgenmax2-tract-8cells_t0000000.vtk \
	axialgenmax2-tract-8cells_t11045160000000.vtk \
	axialgenmax2-tract-8cells_t12623040000000.vtk \
	axialgenmax2-tract-8cells_t14200920000000.vtk \
	axialgenmax2-tract-8cells_t1577880000000.vtk \
	axialgenmax2-tract-8cells_t15778800000000.vtk \
	axialgenmax2-tract-8cells_t3155760000000.vtk \
	axialgenmax2-tract-8cells_t4733640000000.vtk \
	axialgenmax2-tract-8cells_t6311520000000.vtk \
	axialgenmax2-tract-8cells_t7889400000000.vtk \
	axialgenmax2-tract-8cells_t9467280000000.vtk \
	axialgenmax3-tract-8cells_t0000000.vtk \
	axialgenmax3-tract-8cells_t11045160000000.vtk \
	axialgenmax3-tract-8cells_t12623040000000.vtk \
	axialgenmax3-tract-8cells_t14200920000000.vtk \
	axialgenmax3-tract-8cells_t1577880000000.vtk \
	axialgenmax3-tract-8cells_t15778800000000.vtk \
	axialgenmax3-tract-8cells_t3155760000000.vtk \
	axialgenmax3-tract-8cells_t4733640000000.vtk \
	axialgenmax3-tract-8cells_t6311520000000.vtk \
	axialgenmax3-tract-8cells_t7889400000000.vtk \
	axialgenmax3-tract-8cells_t9467280000000.vtk \
	axialmaxwellisotropic_t0.vtk \
	axialmaxwellisotropic_t1.vtk \
	shearelasticisotropic_t0.vtk \
	shearelasticisotropic_t1.vtk \
	shearplanestrain_t0.vtk \
	shearplanestrain_t1.vtk \
	shearmaxwellisotropic_t0.vtk \
	shearmaxwellisotropic_t1.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/3d/hex8
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
