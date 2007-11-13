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
	TestDislocationDyn.py

noinst_DATA = \
	bar.mesh \
	dislocation_dyn.cfg \
	dislocation_slip.spatialdb \
	dislocation_dyn_sliprate.spatialdb \
	dislocation_dyn_sliptime.spatialdb \
	matprops.spatialdb

noinst_TMP = \
	dislocation_dyn_*.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/2d/bar_quad4
export-data: $(noinst_PYTHON) $(noinst_DATA)
	for f in $(noinst_PYTHON) $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_PYTHON) \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
