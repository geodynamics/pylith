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
	box_quad4_100m.exo \
	matprops.spatialdb \
	axialplanestrain.cfg \
	shearplanestrain.cfg \
	dislocation.cfg \
	dislocation_slip.spatialdb \
	dislocation_sliptime.spatialdb

noinst_TMP = \
	axialplanestrain_t0000000.vtk \
	axialplanestrain-statevars-elastic_info.vtk \
	axialplanestrain-statevars-elastic_t0000000.vtk \
	shearplanestrain_t0000000.vtk \
	shearplanestrain-statevars-elastic_info.vtk \
	shearplanestrain-statevars-elastic_t0000000.vtk \
	dislocation_t0000000.vtk \
	dislocation-statevars-elastic_info.vtk \
	dislocation-statevars-elastic_t0000000.vtk \
	dislocation-fault_info.vtk \
	dislocation-fault_t0000000.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/2d/quad4
export-data: $(noinst_PYTHON) $(noinst_DATA)
	for f in $(noinst_PYTHON) $(noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

BUILT_SOURCES = export-data

CLEANFILES = \
	$(export_datadir)/$(noinst_PYTHON) \
	$(export_datadir)/$(noinst_DATA) \
	$(export_datadir)/$(noinst_TMP) \
	*.pyc


# End of file 
