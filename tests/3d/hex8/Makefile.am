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

TESTS = testpylith.py

check_SCRIPTS = testpylith.py

noinst_PYTHON = \
	TestAxialElasticIsotropic.py

noinst_DATA = \
	box.exo \
	matprops.spatialdb \
	matprops-maxwell.spatialdb \
	axialelasticisotropic.cfg \
	axialmaxwellisotropic.cfg \
	axialx.spatialdb \
	axialy.spatialdb \
	axialz.spatialdb \
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
