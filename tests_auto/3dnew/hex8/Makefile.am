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
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

TESTS = testpylith.py

check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestHex8.py \
	TestLgDeformRigidBody.py \
	rigidbody_soln.py \
	rigidbody_gendb.py \
	TestLgDeformTraction.py \
	lgdeformtraction_soln.py \
	testpylith.py

dist_noinst_DATA = \
	geometry.jou \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	lgdeformrigidbody.cfg \
	lgdeformtraction.cfg

noinst_TMP = \
	rigidbody_disp.spatialdb \
	lgdeformrigidbody_t0000000.vtk \
	lgdeformrigidbody-elastic_info.vtk \
	lgdeformrigidbody-elastic_t0000000.vtk \
	lgdeformrigidbody-viscoelastic_info.vtk \
	lgdeformrigidbody-viscoelastic_t0000000.vtk \
	lgdeformtraction_t0000000.vtk \
	lgdeformtraction-elastic_info.vtk \
	lgdeformtraction-elastic_t0000000.vtk \
	lgdeformtraction-viscoelastic_info.vtk \
	lgdeformtraction-viscoelastic_t0000000.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/3dnew/hex8
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc

# End of file 
