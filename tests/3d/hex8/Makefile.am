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
	axialelasticisotropic.cfg \
	axialmaxwellisotropic.cfg \
	axialx.spatialdb \
	axialy.spatialdb \
	axialz.spatialdb \
	matprops.spatialdb \
	matprops-maxwell.spatialdb \
	shearelasticisotropic.cfg \
	shearmaxwellisotropic.cfg \
	shearxyz.spatialdb

noinst_TMP = \
	axialelasticisotropic.vtk \
	axialmaxwellisotropic.vtk \
	shearelasticisotropic.vtk \
	shearmaxwellisotropic.vtk


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
