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
	TestLine2.py \
	TestAxial.py \
	TestDislocation.py

noinst_DATA = \
	bar.mesh \
	matprops.spatialdb \
	axialextension.cfg \
	axialextension_disp.spatialdb \
	dislocation.cfg \
	finalslip.spatialdb \
	sliptime.spatialdb

noinst_TMP = \
	axialextension_t0000000.vtk \
	axialextension-statevars-elastic_info.vtk \
	axialextension-statevars-elastic_t0000000.vtk \
	dislocation-fault_info.vtk \
	dislocation-fault_t0000000.vtk \
	dislocation-statevars-elastic_info.vtk \
	dislocation-statevars-elastic_t0000000.vtk \
	dislocation_t0000000.vtk

TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/1d/line2
export-data: $(noinst_PYTHON) $(noinst_DATA)
	for f in $(noinst_PYTHON) $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_PYTHON) \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP)


# End of file 
