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

dist_noinst_PYTHON = \
	TestQuad4.py \
	TestAxialDisp.py \
	axialdisp_soln.py \
	axialdisp_gendb.py \
	TestShearDisp.py \
	sheardisp_soln.py \
	sheardisp_gendb.py \
	TestDislocation.py \
	dislocation_soln.py \
	dislocation2_soln.py \
	TestLgDeformRigidBody.py \
	rigidbody_soln.py \
	rigidbody_gendb.py \
	TestLgDeformTraction.py \
	lgdeformtraction_soln.py \
	TestFrictionCompression.py \
	friction_compression_soln.py \
	testpylith.py

dist_noinst_DATA = \
	geometry.jou \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	axialdisp.cfg \
	sheardisp.cfg \
	dislocation.cfg \
	dislocation2.cfg \
	lgdeformrigidbody.cfg \
	lgdeformtraction.cfg \
	friction_compression.cfg

noinst_TMP = \
	axial_disp.spatialdb \
	axialdisp_t0000000.vtk \
	axialdisp-elastic_info.vtk \
	axialdisp-elastic_t0000000.vtk \
	shear_disp.spatialdb \
	sheardisp_t0000000.vtk \
	sheardisp-elastic_info.vtk \
	sheardisp-elastic_t0000000.vtk \
	dislocation_t0000000.vtk \
	dislocation-elastic_info.vtk \
	dislocation-elastic_t0000000.vtk \
	dislocation-fault_info.vtk \
	dislocation-fault_t0000000.vtk \
	rigidbody_disp.spatialdb \
	lgdeformrigidbody_t0000000.vtk \
	lgdeformrigidbody-elastic_info.vtk \
	lgdeformrigidbody-elastic_t0000000.vtk \
	lgdeformtraction_t0000000.vtk \
	lgdeformtraction-elastic_info.vtk \
	lgdeformtraction-elastic_t0000000.vtk \
	friction_compression_t0000000.vtk \
	friction_compression-elastic_info.vtk \
	friction_compression-elastic_t0000000.vtk


TESTS_ENVIRONMENT = $(PYTHON)


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests/2d/quad4
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data
clean-local: clean-data
CLEANFILES = *.pyc

# End of file 
