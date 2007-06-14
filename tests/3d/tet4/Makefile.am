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
	TestAxialElasticIsotropic.py \
	TestShearElasticIsotropic.py \
	TestDislocation.py

noinst_DATA = \
	box.exo \
	matprops.spatialdb \
	axialelasticisotropic.cfg \
	axialx.spatialdb \
	axialy.spatialdb \
	axialz.spatialdb \
	shearelasticisotropic.cfg \
	shearxyz.spatialdb \
	twocells.mesh \
	dislocation.cfg \
	dislocation_disp.spatialdb \
	dislocation_sliprate.spatialdb \
	dislocation_slip.spatialdb \
	dislocation_sliptime.spatialdb


noinst_TMP = \
	axialelasticisotropic_t0.vtk \
	axialelasticisotropic_t1.vtk \
	shearelasticisotropic_t0.vtk \
	shearelasticisotropic_t1.vtk \
	dislocation_t0.vtk \
	dislocation_t1.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/3d/tet4
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
