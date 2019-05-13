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
# Copyright (c) 2010-2019 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

subpackage = meshio
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/check.am

SUBDIRS = data

TESTS = testmeshio testvtk

check_PROGRAMS = testmeshio testvtk

# general meshio
testmeshio_SOURCES = \
	test_meshio.cc \
	TestMeshIO.cc \
	TestMeshIOAscii.cc \
	TestMeshIOAscii_Cases.cc \
	TestMeshIOLagrit.cc \
	TestMeshIOLagrit_Cases.cc \
	TestOutputTriggerStep.cc \
	TestOutputTriggerTime.cc

# :TODO: @brad
#TestFieldFilterProject.cc
#TestFieldFilterNone.cc

# VTK data writer
testvtk_SOURCES = \
	test_meshio.cc \
	FieldFactory.cc \
	TestDataWriter.cc \
	TestDataWriterMesh.cc \
	TestDataWriterMaterial.cc \
	TestDataWriterSubmesh.cc \
	TestDataWriterVTK.cc \
	TestDataWriterVTKMesh.cc \
	TestDataWriterVTKMesh_Cases.cc \
	TestDataWriterVTKMaterial.cc \
	TestDataWriterVTKMaterial_Cases.cc \
	TestDataWriterVTKSubmesh.cc \
	TestDataWriterVTKSubmesh_Cases.cc

# :TODO: @brad
# TestDataWriterPoints.cc \
# TestDataWriterVTKPoints.cc \
# TestDataWriterVTKPoints_Cases.cc \
# TestDataWriterFaultMesh.cc \
# TestDataWriterVTKFaultMesh.cc \
# TestDataWriterVTKFaultMesh_Cases.cc \
# TestOutputObserver.cc \
# TestOutputSolnDomain.cc \
# TestOutputSolnBoundary.cc \
# TestOutputSolnPoints.cc 


noinst_HEADERS = \
	TestMeshIO.hh \
	TestMeshIOAscii.hh \
	TestMeshIOLagrit.hh \
	TestOutputTriggerStep.hh \
	TestOutputTriggerTime.hh \
	FieldFactory.hh \
	TestOutputManager.hh \
	TestOutputSolnSubset.hh \
	TestOutputSolnPoints.hh \
	TestVertexFilterVecNorm.hh \
	TestDataWriter.hh \
	TestDataWriterMesh.hh \
	TestDataWriterMaterial.hh \
	TestDataWriterSubmesh.hh \
	TestDataWriterVTK.hh \
	TestDataWriterVTKMesh.hh \
	TestDataWriterVTKMaterial.hh \
	TestDataWriterVTKSubmesh.hh

# :TODO: @brad
# TestDataWriterPoints.hh \
# TestDataWriterVTKPoints.hh
# TestDataWriterFaultMesh.hh \
# TestDataWriterVTKFaultMesh.hh


AM_CPPFLAGS += \
	$(PETSC_SIEVE_FLAGS) $(PETSC_CC_INCLUDES) \
	-I$(PYTHON_INCDIR) $(PYTHON_EGG_CPPFLAGS)

testmeshio_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testmeshio_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

testvtk_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

testvtk_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_CUBIT
  testmeshio_SOURCES += \
	TestExodusII.cc \
	TestMeshIOCubit.cc \
	TestMeshIOCubit_Cases.cc
  noinst_HEADERS += \
	TestExodusII.hh \
	TestMeshIOCubit.hh
  testmeshio_LDADD += -lnetcdf
  testvtk_LDADD += -lnetcdf
endif

if ENABLE_HDF5
  TESTS += testhdf5

  check_PROGRAMS += testhdf5

  # HDF5 data writer
  testhdf5_SOURCES = \
	test_meshio.cc \
	FieldFactory.cc \
	TestDataWriter.cc \
	TestDataWriterMesh.cc \
	TestDataWriterMaterial.cc \
	TestDataWriterSubmesh.cc \
	TestHDF5.cc \
	TestDataWriterHDF5.cc \
	TestDataWriterHDF5Mesh.cc \
	TestDataWriterHDF5Mesh_Cases.cc \
	TestDataWriterHDF5Material.cc \
	TestDataWriterHDF5Material_Cases.cc \
	TestDataWriterHDF5Submesh.cc \
	TestDataWriterHDF5Submesh_Cases.cc \
	TestDataWriterHDF5ExtMesh.cc \
	TestDataWriterHDF5ExtMesh_Cases.cc \
	TestDataWriterHDF5ExtMaterial.cc \
	TestDataWriterHDF5ExtMaterial_Cases.cc \
	TestDataWriterHDF5ExtSubmesh.cc \
	TestDataWriterHDF5ExtSubmesh_Cases.cc


# TestDataWriterHDF5Points.cc \
# TestDataWriterHDF5Point_Cases.cc \
# TestDataWriterHDF5FaultMesh.cc \
# TestDataWriterHDF5FaultMeshCases.cc \
# TestDataWriterHDF5ExtPoints.cc \
# TestDataWriterHDF5ExtPoints_Cases.cc \
# TestDataWriterHDF5ExtFaultMesh.cc \
# TestDataWriterHDF5ExtFaultMeshCases.cc


  noinst_HEADERS += \
	TestHDF5.hh \
	TestDataWriterHDF5.hh \
	TestDataWriterHDF5Mesh.hh \
	TestDataWriterHDF5Material.hh \
	TestDataWriterHDF5Submesh.hh \
	TestDataWriterHDF5ExtMesh.hh \
	TestDataWriterHDF5ExtMaterial.hh \
	TestDataWriterHDF5ExtSubmesh.hh

# TestDataWriterHDF5Points.hh \
# TestDataWriterHDF5FaultMesh.hh \
# TestDataWriterHDF5ExtPoints.hh \
# TestDataWriterHDF5ExtFaultMesh.hh


  testhdf5_LDFLAGS = \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

  testhdf5_LDADD = \
	-lcppunit -ldl \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	-lspatialdata \
	$(PETSC_LIB) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS) \
	-lhdf5

if ENABLE_CUBIT
  testhdf5_LDADD += -lnetcdf
endif

endif


noinst_tmp = \
	mesh1D.txt \
	mesh1Din2D.txt \
	mesh1Din3D.txt \
	mesh2D.txt \
	mesh2Din3D.txt \
	mesh3D.txt


CLEANFILES = $(noinst_tmp)

clean-local: clean-local-tmp
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.vtk *.dat *.dat.info *.h5 *.xmf

# End of file
