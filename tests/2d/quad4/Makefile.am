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
	TestAxialPlaneStrain.py \
	TestShearPlaneStrain.py \
	TestDislocation.py

noinst_DATA = \
	box.exo \
	twocells.mesh \
	matprops.spatialdb \
	axialplanestrain.cfg \
	axialx.spatialdb \
	axialy.spatialdb \
	shearplanestrain.cfg \
	shearplanestrain_free.cfg \
	shearplanestrain_xfixed.cfg \
	shearxy.spatialdb \
	dislocation.cfg \
	dislocation_disp.spatialdb \
	dislocation_sliprate.spatialdb \
	dislocation_slip.spatialdb \
	dislocation_sliptime.spatialdb

noinst_TMP = \
	axialplanestrain_t0.vtk \
	axialplanestrain_t1.vtk \
	shearplanestrain_t0.vtk \
	shearplanestrain_t1.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/2d/quad4
export-data: $(noinst_DATA)
	for f in $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
