# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	apps/__init__.py \
	apps/TestPetscApplication.py \
	apps/TestPyLithApp.py \
	apps/TestEqInfoAPp.py \
	bc/__init__.py \
	bc/TestDirichletTimeDependent.py \
	bc/TestNeumannTimeDependent.py \
	bc/TestAbsorbingDampers.py \
	bc/TestAuxiliarySubfields.py \
	bc/TestZeroDB.py \
	faults/TestFaultCohesive.py \
	faults/TestFaultCohesiveKin.py \
	faults/TestKinSrc.py \
	faults/TestKinSrcConstRate.py \
	faults/TestKinSrcStep.py \
	faults/TestKinSrcRamp.py \
	faults/TestKinSrcBrune.py \
	faults/TestKinSrcLiuCos.py \
	faults/TestKinSrcTimeHistory.py \
	faults/SingleRupture.py \
	materials/__init__.py \
	materials/TestElasticity.py \
	materials/TestHomogeneous.py \
	materials/TestIncompressibleElasticity.py \
	materials/TestMaterial.py \
	materials/TestPoroelasticity.py \
	materials/TestRheologyLinearElasticity.py \
	materials/TestRheologyPoroelasticity.py \
	materials/TestRheologyViscoelasticity.py \
	meshio/__init__.py \
	meshio/TestDataWriter.py \
	meshio/TestDataWriterHDF5.py \
	meshio/TestDataWriterHDF5Ext.py \
	meshio/TestDataWRiterVTK.py \
	meshio/TestXdmf.py


dist_noinst_DATA = \
	meshio/data/tri3.h5 \
	meshio/data/tri3.xmf \
	meshio/data/tri3_vertex.h5 \
	meshio/data/tri3_vertex.xmf \
	meshio/data/tri3_cell.h5 \
	meshio/data/tri3_cell.xmf \
	meshio/data/tri3_points.h5 \
	meshio/data/tri3_points.xmf \
	meshio/data/tri3_points_vertex.h5 \
	meshio/data/tri3_points_vertex.xmf \
	meshio/data/tri3_surf.h5 \
	meshio/data/tri3_surf.xmf \
	meshio/data/tri3_surf_vertex.h5 \
	meshio/data/tri3_surf_vertex.xmf \
	meshio/data/tri3_surf_cell.h5 \
	meshio/data/tri3_surf_cell.xmf \
	meshio/data/quad4.h5 \
	meshio/data/quad4.xmf \
	meshio/data/quad4_vertex.h5 \
	meshio/data/quad4_vertex.xmf \
	meshio/data/quad4_cell.h5 \
	meshio/data/quad4_cell.xmf \
	meshio/data/quad4_points.h5 \
	meshio/data/quad4_points.xmf \
	meshio/data/quad4_points_vertex.h5 \
	meshio/data/quad4_points_vertex.xmf \
	meshio/data/quad4_surf.h5 \
	meshio/data/quad4_surf.xmf \
	meshio/data/quad4_surf_vertex.h5 \
	meshio/data/quad4_surf_vertex.xmf \
	meshio/data/quad4_surf_cell.h5 \
	meshio/data/quad4_surf_cell.xmf \
	meshio/data/tet4.h5 \
	meshio/data/tet4.xmf \
	meshio/data/tet4_vertex.h5 \
	meshio/data/tet4_vertex.xmf \
	meshio/data/tet4_cell.h5 \
	meshio/data/tet4_cell.xmf \
	meshio/data/tet4_points.h5 \
	meshio/data/tet4_points.xmf \
	meshio/data/tet4_points_vertex.h5 \
	meshio/data/tet4_points_vertex.xmf \
	meshio/data/tet4_surf.h5 \
	meshio/data/tet4_surf.xmf \
	meshio/data/tet4_surf_vertex.h5 \
	meshio/data/tet4_surf_vertex.xmf \
	meshio/data/tet4_surf_cell.h5 \
	meshio/data/tet4_surf_cell.xmf \
	meshio/data/hex8.h5 \
	meshio/data/hex8.xmf \
	meshio/data/hex8_vertex.h5 \
	meshio/data/hex8_vertex.xmf \
	meshio/data/hex8_cell.h5 \
	meshio/data/hex8_cell.xmf \
	meshio/data/hex8_points.h5 \
	meshio/data/hex8_points.xmf \
	meshio/data/hex8_points_vertex.h5 \
	meshio/data/hex8_points_vertex.xmf \
	meshio/data/hex8_surf.h5 \
	meshio/data/hex8_surf.xmf \
	meshio/data/hex8_surf_vertex.h5 \
	meshio/data/hex8_surf_vertex.xmf \
	meshio/data/hex8_surf_cell.h5 \
	meshio/data/hex8_surf_cell.xmf

noinst_TMP = \
	meshio/tri3.xmf \
	meshio/tri3_vertex.xmf \
	meshio/tri3_cell.xmf \
	meshio/tri3_points.xmf \
	meshio/tri3_points_vertex.xmf \
	meshio/tri3_surf.xmf \
	meshio/tri3_surf_vertex.xmf \
	meshio/tri3_surf_cell.xmf \
	meshio/quad4.xmf \
	meshio/quad4_vertex.xmf \
	meshio/quad4_cell.xmf \
	meshio/quad4_points.xmf \
	meshio/quad4_points_vertex.xmf \
	meshio/quad4_surf.xmf \
	meshio/quad4_surf_vertex.xmf \
	meshio/quad4_surf_cell.xmf \
	meshio/tet4.xmf \
	meshio/tet4_vertex.xmf \
	meshio/tet4_cell.xmf \
	meshio/tet4_points.xmf \
	meshio/tet4_points_vertex.xmf \
	meshio/tet4_surf.xmf \
	meshio/tet4_surf_vertex.xmf \
	meshio/tet4_surf_cell.xmf \
	meshio/hex8.xmf \
	meshio/hex8_vertex.xmf \
	meshio/hex8_cell.xmf \
	meshio/hex8_points.xmf \
	meshio/hex8_points_vertex.xmf \
	meshio/hex8_surf.xmf \
	meshio/hex8_surf_vertex.xmf \
	meshio/hex8_surf_cell.xmf

# 'export' the input files by performing a mock install
export_datadir = $(abs_top_builddir)/tests/pytests
export-data: $(dist_noinst_DATA)
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir)/$$f; done; fi

clean-data:
	if [ "X$(abs_top_srcdir)" != "X$(abs_top_builddir)" ]; then for f in $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi

BUILT_SOURCES = export-data
clean-local: clean-data


# End of file 
